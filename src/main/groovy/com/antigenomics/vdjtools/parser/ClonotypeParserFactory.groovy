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

import com.antigenomics.vdjtools.SampleJ
import com.antigenomics.vdjtools.Software

class ClonotypeParserFactory {
    public static ClonotypeParser getParser(Software software, SampleJ sample) {
        switch (software) {
            case Software.MiTcr:
                return new MiTcrParser(sample)
            case Software.IgBlast:
                return new IgBlastParser(sample)
            case Software.Simple:
                return new SimpleParser(sample)
            case Software.MiGec:
                return new MiGecParser(sample)
            default:
                throw new UnsupportedOperationException("Don't know how to parse $software data")
        }
    }
}
