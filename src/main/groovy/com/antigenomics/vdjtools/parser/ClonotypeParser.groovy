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

package com.antigenomics.vdjtools.parser

import com.antigenomics.vdjtools.ClonotypeJ
import com.antigenomics.vdjtools.SampleJ
import com.antigenomics.vdjtools.Software

public abstract class ClonotypeParser implements Iterable<ClonotypeJ> {
    public static ClonotypeParser create(InputStream inputStream, Software software, SampleJ sample) {
        ClonotypeParser parser
        def reader = new BufferedReader(new InputStreamReader(inputStream))
        def innerIter = reader.iterator()

        switch (software) {
            case Software.MiTcr:
                parser = new MiTcrParser(innerIter: innerIter, software: software, sample: sample)
                break
            case Software.IgBlast:
                parser = new IgBlastParser(innerIter: innerIter, software: software, sample: sample)
                break
            case Software.Simple:
                parser = new SimpleParser(innerIter: innerIter, software: software, sample: sample)
                break
            case Software.MiGec:
                parser = new MiGecParser(innerIter: innerIter, software: software, sample: sample)
                break
            default:
                throw new UnsupportedOperationException("Don't know how to parse $software data")
        }

        (1..software.headerLineCount).each { innerIter.next() }

        return parser
    }

    protected final Software software
    protected final Iterator<String> innerIter
    protected final SampleJ sample

    protected abstract ClonotypeJ parse(String clonotypeString)

    /**
     * Internal util, extracts most probable allele
     */
    protected static List<String> extractVDJ(List<String> vdj) {
        vdj.collect {
            def major = it.split(",")[0]
            major = major.split("\\*")[0] // trim allele if present
            major = major.replaceAll("\"", "")  // zap characters introduced by opening file in Excel
            major.length() > 0 ? major : "."
        }
    }

    @Override
    Iterator<ClonotypeJ> iterator() {
        [hasNext: {
            innerIter.hasNext()
        }, next : {
            parse(innerIter.next())
        }] as Iterator
    }
}
