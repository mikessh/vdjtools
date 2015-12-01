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

package com.antigenomics.vdjtools.join;

import java.util.Collections;
import java.util.Set;

/**
 * A filter that retains or filters out all clonotypes that were detected in a specified sample(s).
 * Decision is made according to the number of times the clonotype was detected in the sample set,
 * not accounting for convergent variant count.
 */
public class SampleSpecificOccurrenceJoinFilter implements JoinFilter {
    private final int occurrenceThreshold;
    private final Set<String> sampleIds;
    private final boolean enrichment;

    /**
     * Creates a sample-specific filter. Clonotype occurrences are counted in a specified sample set.
     * Decision to retain or filter joint clonotype is made according to occurrence count.
     *
     * @param sampleIds           identifiers (sample.id from metadata) of samples.
     * @param occurrenceThreshold inclusive lower bound on number of times the clonotype should be detected in specified sample set to be marked as "detected".
     * @param enrichment          if set to true will retain clonotypes that were detected, will filter them out otherwise.
     */
    public SampleSpecificOccurrenceJoinFilter(Set<String> sampleIds,
                                              int occurrenceThreshold, boolean enrichment) {
        this.occurrenceThreshold = occurrenceThreshold;
        this.sampleIds = sampleIds;
        this.enrichment = enrichment;
    }

    public int getOccurrenceThreshold() {
        return occurrenceThreshold;
    }

    public Set<String> getSampleIds() {
        return Collections.unmodifiableSet(sampleIds);
    }

    public boolean isEnrichment() {
        return enrichment;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean pass(JointClonotype jointClonotype) {
        int detectionCounter = 0;
        for (int i = 0; i < jointClonotype.getParent().getNumberOfSamples(); i++) {
            if (sampleIds.contains(jointClonotype.getParent().getSample(i).getSampleMetadata().getSampleId()) &&
                    jointClonotype.present(i) && ++detectionCounter == occurrenceThreshold) return enrichment;
        }
        return !enrichment;
    }
}
