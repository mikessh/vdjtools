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

import com.antigenomics.vdjtools.ClonotypeContainer
import com.antigenomics.vdjtools.util.CommonUtil

class CdrPwmGrid {
    private final Map<Cell, CdrPwm> pwmGrid = Collections.synchronizedMap(new HashMap<>())

    public void update(ClonotypeContainer clonotypes) {
        clonotypes.each {
            def pwm = getAtOrCreate(it.v, it.cdr3aa.length())
            pwm.update(it)
        }
    }

    private CdrPwm getAtOrCreate(String v, int length) {
        def cell = new Cell(v, length)
        def pwm = pwmGrid[cell]
        if (!pwm)
            pwmGrid.put(cell, pwm = new CdrPwm(length))
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

    public Collection<Cell> filterCells(int minCount,
                                        double freqThreshold) {
        pwmGrid.findAll { it.value.count >= minCount && it.value.totalFreq >= freqThreshold }.keySet()
    }

    public Map<Cell, CdrPattern> compilePatterns(int minCount,
                                                 double freqThreshold,
                                                 double positionalRatioThreshold) {
        filterCells(minCount, freqThreshold).collectEntries {
            [(it): pwmGrid[it].extractConsensus(positionalRatioThreshold)]
        }
    }

    class Cell {
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


    public static final String HEADER = "#v\tlen\tpos\t" + CommonUtil.AAS.collect().join("\t")

    public String toString(int minCount, double freqThreshold) {
        filterCells(minCount, freqThreshold).collect { Cell cell ->
            def pwm = pwmGrid[cell]
            pwm.collect { CdrPwm.Row row ->
                [cell.v, cell.length, row.pos, row.freqsNorm.collect()].flatten().join("\t")
            }.join("\n")
        }.join("\n")
    }

    @Override
    public String toString() {
        toString(0, 0)
    }
}
