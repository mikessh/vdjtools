/*
 * Copyright 2013-2014 Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Last modified on 23.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.pwm

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.ClonotypeContainer
import com.antigenomics.vdjtools.util.CommonUtil
import com.antigenomics.vdjtools.util.ExecUtil
import groovyx.gpars.GParsPool

class CdrPwmGrid {
    private final Map<Cell, CdrPwm> pwmGrid = Collections.synchronizedMap(new HashMap<>())

    public static CdrPwmGrid CONTROL = fromInputStream(CommonUtil.resourceStreamReader("pwm/healthy70.txt"))

    public static CdrPwmGrid fromInputStream(InputStreamReader reader) {
        def cdrPwmGrid = new CdrPwmGrid()
        def firstLine
        while ((firstLine = reader.readLine())) {
            if (!firstLine.startsWith("#")) {
                def splitLine = firstLine.split("\t")
                def v = splitLine[0], length = splitLine[1].toInteger(),
                    uniq = splitLine[2].toInteger(), freq = splitLine[3].toDouble()
                def cell = new Cell(v, length)
                def pwm = new double[length][CommonUtil.AAS.length]
                for (int i = 0; i < length; i++) {
                    splitLine = reader.readLine().split("\t")
                    splitLine[3..-1].eachWithIndex { it, j -> // exclude v/len/pos cols
                        pwm[i][j] = it.toDouble()
                    }
                }
                def cdrPwm = new CdrPwm(cell, pwm, uniq, freq)
                cdrPwmGrid.pwmGrid.put(cell, cdrPwm)
            }
        }
        cdrPwmGrid
    }

    public void update(ClonotypeContainer clonotypes) {
        GParsPool.withPool ExecUtil.THREADS, {
            clonotypes.eachParallel { Clonotype clonotype ->
                if (clonotype.coding) {
                    def pwm = getAtOrCreate(clonotype.v, clonotype.cdr3aa.length())
                    pwm.update(clonotype)
                }
            }
        }
    }

    private CdrPwm getAtOrCreate(String v, int length) {
        def cell = new Cell(v, length)
        def pwm = pwmGrid[cell]
        if (!pwm)
            pwmGrid.put(cell, pwm = new CdrPwm(cell))
        pwm
    }

    public CdrPwm getAt(Cell cell) {
        pwmGrid[cell]
    }

    public CdrPwm getAt(String v, int length) {
        getAt(new Cell(v, length))
    }

    public Collection<CdrPwm> getAt(String v) {
        pwmGrid.findAll { it.key.v == v }.values()
    }

    public Collection<Cell> collectCells() {
        Collections.unmodifiableCollection(pwmGrid.keySet())
    }

    public Set<String> getVSegments() {
        new HashSet<String>(collectCells().collect { it.v })
    }

    public Set<Integer> getLengths() {
        new HashSet<Integer>(collectCells().collect { it.length })
    }

    public Collection<Cell> filterCells(int minCount,
                                        double freqThreshold) {
        pwmGrid.findAll { it.value.div >= minCount && it.value.freq >= freqThreshold }.keySet()
    }

    public Map<Cell, CdrPattern> compilePatterns(int minCount,
                                                 double freqThreshold,
                                                 double positionalRatioThreshold) {
        filterCells(minCount, freqThreshold).collectEntries {
            [(it): pwmGrid[it].extractConsensus(positionalRatioThreshold)]
        }
    }

    static class Cell {
        public final int length
        public final String v

        Cell(String v, int length) {
            this.length = length
            this.v = v
        }

        @Override
        boolean equals(o) {
            Cell that = (Cell) o

            length == that.length && v == that.v
        }

        @Override
        int hashCode() {
            31 * length + v.hashCode()
        }
    }

    public static final String HEADER = "#v\tlen\tpos\tuniq\tfreq\t" + CommonUtil.AAS.collect().join("\t")

    public String toString(int minCount, double freqThreshold, boolean normalize, boolean correct) {
        filterCells(minCount, freqThreshold).collect { Cell cell ->
            def pwm = pwmGrid[cell]
            (0..<cell.length).collect { int pos ->
                [cell.v, cell.length, pos,
                 pwm.div, pwm.freq,
                 pwm.getNormalizedFreqs(pos, normalize ? CONTROL : null, correct)].flatten().join("\t")
            }.join("\n")
        }.join("\n")
    }

    @Override
    public String toString() {
        toString(1, 0.01, true, false)
    }
}
