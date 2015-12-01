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


package com.antigenomics.vdjtools.io.parser

import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.misc.CommonUtil

import static com.antigenomics.vdjtools.misc.CommonUtil.toUnifiedCdr3Aa

/**
 * A clonotype parser implementation that handles output from IgBlastWrapper software.
 * {@url https://github.com/mikessh/igblastwrp}
 */
public class MigMapParser extends ClonotypeStreamParser {
    /**
     * {@inheritDoc}
     */
    protected MigMapParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.MigMap, sample)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        /*
             0     1      2  3  4  5       6           7    8     9    10    11   12    13
            "freq\tcount\tv\td\tj\tcdr3nt\tcdr3aa\t" + FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4 +
               14               15                  16    17      18    19 
            "\tcdr.insert.qual\tmutations.qual\t" + vEnd, dStart, dEnd, jStart + "\t" +
            20    21     22     23
            vDel, dDel5, dDel3, jDel +
            24      25       26       27
            "pol.v\tpol.d.5\tpol.d.3\tpol.j" +
               28        29        30       31        32
            "\thas.cdr3\tin.frame\tno.stop\tcomplete\tcanonical"
         */

        def splitString = clonotypeString.split("\t")

        def freq = splitString[0].toDouble()
        def count = splitString[1].toInteger()
        def cdr3nt = splitString[5] == "." ? "" : splitString[5],
            cdr3aa = splitString[6] == "." ? "" : toUnifiedCdr3Aa(splitString[6])

        String v, d, j
        (v, d, j) = CommonUtil.extractVDJ(splitString[2..4])

        boolean inFrame, noStop, isComplete
        (inFrame, noStop, isComplete) = splitString[29..31].collect { it.toBoolean() }

        def segmPoints = splitString[16..19].collect { it.toInteger() } as int[]

        def clonotype = new Clonotype(sample,
                count, freq,
                segmPoints, v, d, j,
                cdr3nt, cdr3aa,
                inFrame, noStop, isComplete)

        return clonotype
    }
}
