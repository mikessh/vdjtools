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

package com.antigenomics.vdjtools.clustering.spectratype

import com.antigenomics.vdjtools.clustering.ClonotypeCluster
import com.antigenomics.vdjtools.clustering.SingleSampleClustering
import com.antigenomics.vdjtools.sample.Sample

class Spectratype implements SingleSampleClustering {
    private final Map<String, SpectratypePeak> peaks = new HashMap<>()
    final Sample parentSample
    int clones = 0
    double freq = 0

    Spectratype(Sample parentSample) {
        this.parentSample = parentSample

        parentSample.each { clonotype ->
            def tempPeak = new SpectratypePeak(clonotype)

            def peak = peaks[tempPeak.signature]
            if (!peak) {
                peak = tempPeak
                peaks.put(peak.signature, peak)
            }

            peak.append(clonotype)

            clones++
            freq += clonotype.freq
        }
    }

    //private List<SpectratypePeak> sortedPeaks = null

    //List<SpectratypePeak> getSortedPeaks() {
    //    sortedPeaks ?: (sortedPeaks = peaks.sort { it.key }.collect { it.value })
    //}

    //int getNumberOfPeaks() {
    //    peaks.size()
    //}

    @Override
    String toString() {
        "#" + SpectratypePeak.HEADER + "\n" + clusters.join("\n")
    }

    @Override
    List<ClonotypeCluster> getClusters() {
        peaks.values()
    }

    @Override
    int getNumberOfClusters() {
        peaks.size()
    }
}
