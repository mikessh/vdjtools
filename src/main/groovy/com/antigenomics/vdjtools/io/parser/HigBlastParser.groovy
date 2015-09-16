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

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.CommonUtil

import static com.antigenomics.vdjtools.util.CommonUtil.toUnifiedCdr3Aa

/**
 * A clonotype parser implementation that handles output from IgBlastWrapper software.
 * {@url https://github.com/mikessh/igblastwrp}
 */
public class HigBlastParser extends ClonotypeStreamParser {
    /**
     * {@inheritDoc}
     */
    protected HigBlastParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.HigBlast, sample)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        /*
        0	1	2	3
        #reads_count	reads_percent	events_count	events_percent
        4	5	6	7	8	9
        cdr1nt	cdr2nt	cdr3nt	cdr1aa	cdr2aa	cdr3aa
        10	11	12
        inFrame	noStop	complete
        13	14	15
        vSegment	dSegment	jSegment
        16	17	18
        cdr1q	cdr2q	cdr3q
        19
        mutations
         */

        def splitString = clonotypeString.split("\t")

        def count = splitString[1].toInteger()
        def freq = splitString[0].toDouble()
        def cdr3nt = splitString[5] == "." ? "" : splitString[5],
            cdr3aa = splitString[6] == "." ? "" : toUnifiedCdr3Aa(splitString[6])

        String v, d, j
        (v, d, j) = CommonUtil.extractVDJ(splitString[2..4])

        boolean inFrame, noStop, isComplete
        (inFrame, noStop, isComplete) = splitString[25..27].collect { it.toBoolean() }

        def segmPoints = splitString[16..19].collect { it.toInteger() } as int[]

        def clonotype = new Clonotype(sample,
                count, freq,
                segmPoints, v, d, j,
                cdr3nt, cdr3aa,
                inFrame, noStop, isComplete)

        //if (splitString[19] != ".")
        //    mutations.addAll(splitString[19].split("\\|").collect { String mutString ->
        //        parseIgBlastMutation(mutString, clonotype)
        //    })

        return clonotype
    }

    /*
    private static Mutation parseIgBlastMutation(String mutationString, Clonotype parent) {
        def splitString = mutationString.split(",")

        def ntString = splitString[1]
        def splitNTString = ntString.split("[:>]")
        def ntPos = splitNTString[0].toInteger()
        char fromNt = splitNTString[1]
        char toNt = splitNTString[2]

        def aaString = splitString[2]
        def splitAAString = aaString.split("[:>]")
        def aaPos = splitAAString[0].toInteger()
        char fromAa = splitAAString[1]
        char toAa = splitAAString[2]

        def region = splitString[3]

        new Mutation(region,
                ntPos, aaPos,
                fromAa, toAa, fromNt, toNt,
                true, parent)
    }*/
}
