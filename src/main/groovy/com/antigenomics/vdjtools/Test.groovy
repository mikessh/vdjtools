package com.antigenomics.vdjtools

import com.antigenomics.vdjtools.table.HypermutationsByRegion

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

println HypermutationsByRegion.HEADER


def a = ["Reads", "Clones"], b = ["Silent", "Replacement"], c = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "N/A"]

def str = new String[2][2][6]

(0..1).each { i ->
    (0..1).each { j ->
        (0..5).each { k ->
            str[i][j][k] = a[i] + b[j] + c[k]
        }
    }
}

println str.collect().flatten().join("\t")