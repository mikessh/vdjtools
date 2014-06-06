package com.antigenomics.vdjtools.imgt

import java.util.regex.Matcher
import java.util.regex.Pattern

/**
 Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */
class ImgtToMigecParser {
    final static int IMGT_V_REF = 312
    final static String CysRegex = /TG[TC]/, PheTrpRegex = /(?:TGG|TT[TC])(?:GG[ATGC]|GC[ATGC])...GG[ATGC]/
    final boolean onlyFunctional, onlyMajorAllele

    final List<ImgtRecord> failedVReferencePoint = new ArrayList<>(),
                           failedJReferencePoint = new ArrayList<>(),
                           otherSegment = new ArrayList<>()

    ImgtToMigecParser(boolean onlyFunctional, boolean onlyMajorAllele) {
        this.onlyFunctional = onlyFunctional
        this.onlyMajorAllele = onlyMajorAllele
    }

    static String getGene(String id) {
        id.substring(0, 4)
    }

    static String sequenceNoGaps(ImgtRecord record) {
        record.sequence.replaceAll("\\.", "")
    }

    static boolean majorAllele(ImgtRecord record) {
        record.allele == "01"
    }

    static boolean functional(ImgtRecord record) {
        record.functionality == "F"
    }

    static int getVReferencePoint(ImgtRecord imgtRecord) {
        String sequenceWithGaps = imgtRecord.sequence
        if (IMGT_V_REF > sequenceWithGaps.length())
            return -1
        String codon = sequenceWithGaps.substring(IMGT_V_REF - 3, IMGT_V_REF)
        if (!(codon =~ CysRegex))
            return -1
        int imgtGapsCount = sequenceWithGaps.substring(0, IMGT_V_REF).count(".")
        IMGT_V_REF - imgtGapsCount
    }

    static int getJReferencePoint(ImgtRecord imgtRecord) {
        String sequenceWithoutGaps = sequenceNoGaps(imgtRecord) // just in case
        Matcher matcher = Pattern.compile(PheTrpRegex).matcher(sequenceWithoutGaps);
        matcher.find() ? matcher.start() - 1 : -1
    }

    MigecSegmentRecord createRecord(ImgtRecord imgtRecord, String segment, int referencePoint) {
        new MigecSegmentRecord(imgtRecord.species, getGene(imgtRecord.fullId),
                segment, imgtRecord.fullId, sequenceNoGaps(imgtRecord),
                referencePoint)
    }

    MigecSegmentRecord parseRecord(ImgtRecord imgtRecord) {
        if ((!onlyFunctional || functional(imgtRecord)) &&
                (!onlyMajorAllele || majorAllele(imgtRecord))) {
            switch (imgtRecord.type) {
                case "V-REGION":
                    int vReferencePoint = getVReferencePoint(imgtRecord)
                    if (vReferencePoint < 0) {
                        failedVReferencePoint.add(imgtRecord)
                        return null
                    }
                    return createRecord(imgtRecord, "Variable", vReferencePoint)
                case "D-REGION":
                    return createRecord(imgtRecord, "Diversity", -1)
                case "J-REGION":
                    int jReferencePoint = getJReferencePoint(imgtRecord)
                    if (jReferencePoint < 0) {
                        failedJReferencePoint.add(imgtRecord)
                        return null
                    }
                    return createRecord(imgtRecord, "Joining", jReferencePoint)
            }
        }
        otherSegment.add(imgtRecord)
        null
    }
}
