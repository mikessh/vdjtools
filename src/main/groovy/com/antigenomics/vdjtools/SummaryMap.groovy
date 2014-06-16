package com.antigenomics.vdjtools

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
class SummaryMap {
    final Map<String, double[]> map = new HashMap<>()
    final int nCounters
    final String header

    SummaryMap(int nCounters, String header) {
        this.header = header
        this.nCounters = nCounters
    }

    void put(String key, List counters) {
        double[] oldCounters = map[key]
        if (oldCounters == null)
            map.put(key, oldCounters = new double[nCounters])

        for (int i = 0; i < nCounters; i++)
            oldCounters[i] += counters[i]
    }

    void writeToFile(String fileName) {
        new File(fileName).withPrintWriter { pw ->
            pw.println(header)
            map.sort { -it.value[0] }.each {
                pw.println(it.value.collect().join("\t") + "\t" + it.key)
            }
        }
    }
}
