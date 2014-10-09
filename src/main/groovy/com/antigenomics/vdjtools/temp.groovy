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


package com.antigenomics.vdjtools

def samples1 = [], samples2 = []

new File(args[0]).splitEachLine("\t") {
    if (it[1] == "1")
        samples1.add(it[0])
    else
        samples1.add(it[1])
}

def aa1 = new HashMap(), aa2 = new HashMap()

samples1.each { f ->
    new File(f).splitEachLine("\t") {
        aa1.put(it[5], (aa1[it[5]] ?: 0) + 1)
    }
}


samples1.each { f ->
    new File("c1_${f}").withPrintWriter { pw ->
        new File(f).splitEachLine("\t") {
            def count = aa1[it[5]]
            if (count && count > 1)
                pw.println(it.collect().join("\t"))
        }
    }
}

samples2.each { f ->
    new File(f).splitEachLine("\t") {
        aa2.put(it[5], (aa2[it[5]] ?: 0) + 1)
    }
}


samples2.each { f ->
    new File("c2_${f}").withPrintWriter { pw ->
        new File(f).splitEachLine("\t") {
            def count = aa2[it[5]]
            if (count && count > 1)
                pw.println(it.collect().join("\t"))
        }
    }
}