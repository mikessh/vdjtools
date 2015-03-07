/*
 * Copyright 2013-2015 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 */


package com.antigenomics.vdjtools.io.parser

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.CommonUtil

/**
 * A clonotype parser implementation that handles Adaptive Biotechnologies (tm) immunoSEQ (tm) assay
 * output format, see
 * {@url http://www.adaptivebiotech.com/content/immunoseq-0}
 */
class AdaptiveParser extends ClonotypeStreamParser {
    /**
     * {@inheritDoc}
     */
    public AdaptiveParser(Iterator<String> innerIter, Software software, Sample sample) {
        super(innerIter, software, sample)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        /*
             $id   | content                    | description
            -------|----------------------------|------------------------------------------------------------------------
             00    | nucleotide sequence        | this is raw sequence of a read. It could contain only 5' part of CDR3
             01    | CDR3 amino acid sequence   | standard AA sequence or blank, in case read was too short to cover CDR3
             02    | Count                      | default
             03    | Fraction                   | fraction that also accounts for incomplete reads
             04-05 | unused                     |
             06    | V segments                 | V "family", non-conventional naming
             07-12 | unused                     |
             13    | D segments                 | D "family", non-conventional naming
             14-19 | unused                     |
             20    | J segments                 | J "family", non-conventional naming
             21-31 | unused                     |
             32    | CDR3 start                 | used to extract CDR3 nucleotide sequence from $00
             33    | V end                      |
             34    | D start                    |
             35    | D end                      |
             36    | J start                    |
             37    | unused                     |
             38    | status                     | In/Out/Stop
             39+   | unused                     |
             
             - note that there is no "J end" (perhaps due to J segment identification issues), so
              here we use $32 + length($01) * 3
          */

        def splitString = clonotypeString.split("\t")

        def count = splitString[2].toInteger()
        def freq = splitString[3].toDouble()

        def cdr3start = splitString[32].toInteger()

        String cdr1nt = null, cdr2nt = null, cdr3nt, cdr1aa = null, cdr2aa = null, cdr3aa
        cdr3aa = splitString[1]
        cdr3nt = cdr3aa.length() > 0 ? splitString[0].substring(cdr3start, cdr3start + cdr3aa.length() * 3) : ""


        String v, d, j
        (v, d, j) = CommonUtil.extractVDJ(splitString[[6, 13, 20]])

        def status = splitString[38].toUpperCase()

        boolean inFrame = status != "OUT", noStop = status != "STOP", isComplete = cdr3aa.length() > 0

        def segmPoints = [
                splitString[33].toInteger() - 1,
                splitString[34].toInteger(),
                splitString[35].toInteger() - 1,
                splitString[36].toInteger()] as int[]

        new Clonotype(sample, count, freq,
                segmPoints, v, d, j,
                cdr1nt, cdr2nt, cdr3nt,
                cdr1aa, cdr2aa, cdr3aa,
                inFrame, noStop, isComplete,
                new HashSet<>())
    }
}
