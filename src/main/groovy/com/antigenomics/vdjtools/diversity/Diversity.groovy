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

package com.antigenomics.vdjtools.diversity

class Diversity {
    public final long mean, std, n
    public final boolean resampled
    public final DiversityType type
    public final String estimatorName

    public static Diversity DUMMY = new Diversity(-1, -1, -1, DiversityType.Unknown, false, "DUMMY")

    Diversity(double mean, double std, long n,
              DiversityType type, boolean resampled, String estimatorName) {
        this((long) mean, (long) std, n, type, resampled, estimatorName)
    }

    Diversity(long mean, long std, long n,
              DiversityType type, boolean resampled, String estimatorName) {
        this.mean = mean
        this.std = std
        this.n = n
        this.resampled = resampled
        this.type = type
        this.estimatorName = estimatorName
        this.HEADER = [estimatorName + "_mean", estimatorName + "_std"]
    }

    public final String HEADER

    @Override
    public String toString() {
        [mean, std].join("\t")
    }
}
