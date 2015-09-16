/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */

package com.antigenomics.vdjtools.profile

import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.ClonotypeContainer
import com.antigenomics.vdjtools.util.CommonUtil
import com.antigenomics.vdjtools.util.ExecUtil
import groovyx.gpars.GParsPool

class CdrPwmGrid {
    private double freq = 0
    private int div = 0
    private final Map<Cell, CdrPwm> pwmGrid = Collections.synchronizedMap(new HashMap<>())

    public static CdrPwmGrid CONTROL = fromInputStream(CommonUtil.resourceStreamReader("profile/healthy70.txt"))

    private static double convertFreq(String s) {
        double x = s.toDouble()

        if (x > 1.001 || x < -0.001) {
            println "[ERROR] It is highly likely that you're trying to load a PWM grid " +
                    "that was created without --raw parameter. In such case entropy-scaled " +
                    "values are saved, which would not be properly handled by current routines. " +
                    "Aborting..."
            System.exit(-1)
        }

        x
    }

    public static CdrPwmGrid fromInputStream(InputStreamReader reader) {
        def cdrPwmGrid = new CdrPwmGrid()
        def firstLine
        while ((firstLine = reader.readLine())) {
            if (!firstLine.startsWith("#")) {
                // #v	len	pos	div	freq
                def splitLine = firstLine.split("\t")
                def v = splitLine[0], length = splitLine[1].toInteger(),
                    div = splitLine[3].toInteger(), freq = splitLine[4].toDouble()
                def cell = new Cell(v, length)
                def pwm = new double[length][CommonUtil.AAS.length]

                cdrPwmGrid.freq += freq
                cdrPwmGrid.div += div

                // first line
                splitLine[5..-1].eachWithIndex { String it, int j -> // exclude v/len/pos cols
                    pwm[0][j] = convertFreq(it) * freq
                }

                for (int i = 1; i < length; i++) {
                    splitLine = reader.readLine().split("\t")
                    splitLine[5..-1].eachWithIndex { String it, int j -> // exclude v/len/pos cols
                        pwm[i][j] = convertFreq(it) * freq
                    }
                }

                def cdrPwm = new CdrPwm(cdrPwmGrid, cell, pwm, div, freq)
                cdrPwmGrid.pwmGrid.put(cell, cdrPwm)
            }
        }
        cdrPwmGrid
    }

    public void update(ClonotypeContainer clonotypes) {
        this.freq += clonotypes.freq
        this.div += clonotypes.diversity
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
            pwmGrid.put(cell, pwm = new CdrPwm(this, cell))
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

    double getFreq() {
        freq
    }

    int getDiv() {
        div
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

    public static final String HEADER = "#v\tlen\tpos\tdiv\tfreq\t" + CommonUtil.AAS.collect().join("\t")

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

    public String toStringRaw() {
        pwmGrid.collect {
            def cell = it.key
            def pwm = it.value
            (0..<cell.length).collect { int pos ->
                [cell.v, cell.length, pos,
                 pwm.div, pwm.freq,
                 pwm.getFreqs(pos)].flatten().join("\t")
            }.join("\n")
        }.join("\n")
    }

    @Override
    public String toString() {
        toString(1, 0.01, true, false)
    }
}
