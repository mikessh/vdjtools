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

package com.antigenomics.vdjtools.sample

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.metadata.SampleMetadata

class Sample implements Iterable<Clonotype> {
    final SampleMetadata sampleMetadata
    final List<Clonotype> clonotypes
    private Long count = null
    private Integer diversityCDR3NT = null, diversityCDR3AA = null
    private Double freq = null

    Sample(SampleMetadata sampleMetadata, List<Clonotype> clonotypes) {
        this.sampleMetadata = sampleMetadata
        this.clonotypes = clonotypes
    }

    Sample(String sampleId, List<Clonotype> clonotypes) {
        this.sampleMetadata = SampleMetadata.create(sampleId)
        this.clonotypes = clonotypes
    }

    Sample top(int numberOfClonotypes) {
        numberOfClonotypes = Math.min(numberOfClonotypes + 1, clonotypes.size())
        new Sample(sampleMetadata, clonotypes.sort { -it.freq }.subList(0, numberOfClonotypes - 1))
    }

    Sample functional() {
        def sample = new Sample(sampleMetadata, clonotypes.findAll {
            it.inFrame && it.noStop
        })
        sample.renormalize()
        sample
    }

    Sample clone() {
        new Sample(sampleMetadata, clonotypes)
    }

    int getDiv() {
        clonotypes.size()
    }

    int getDiversityCDR3NT() {
        diversityCDR3NT ?: (diversityCDR3NT = new HashSet<String>(clonotypes.collect { it.cdr3nt }).size())
    }

    int getDiversityCDR3AA() {
        diversityCDR3AA ?: (diversityCDR3AA = new HashSet<String>(clonotypes.collect { it.cdr3aa }).size())
    }

    long getCount() {
        count ?: (count = (long) clonotypes.sum { it.count })
    }

    double getFreq() {
        freq ?: (freq = (double) clonotypes.sum { it.freq })
    }

    void renormalize() {
        clonotypes.each {
            it.freq = it.count / (double) getCount()
        }
    }

    //@Override
    //String toString() {
    //    metadata.toString()
    //}


    void sort() {
        clonotypes.sort { -it.freq }
    }

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        Sample sample = (Sample) o

        if (sampleMetadata != sample.sampleMetadata) return false

        return true
    }

    @Override
    int hashCode() {
        return sampleMetadata.hashCode()
    }

    @Override
    Iterator<Clonotype> iterator() {
        clonotypes.iterator()
    }

    //
    // Table output
    //

    void print(PrintWriter printWriter, Software software, boolean addHeader) {
        if (addHeader)
            printWriter.println(software.header)
        clonotypes.each {
            it.print(printWriter, software)
        }
    }
}
