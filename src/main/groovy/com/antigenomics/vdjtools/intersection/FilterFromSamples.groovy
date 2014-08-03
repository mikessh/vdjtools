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

package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.Software

def cli = new CliBuilder(usage: "FilterFromSamples [options] filter_sample sample1 [sample2 ...] output_prefix")
cli.h("display help message")
cli.i(longOpt: "intersect-type", argName: "string1,string2,..", args: 1,
        "Comma-separated list of intersection types to apply. " +
                "Allowed values: $IntersectionType.allowedNames. " +
                "Will use '$IntersectionType.AminoAcid.shortName' by default.")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.n(longOpt: "negative", "Will report clonotypes present in filter_sample. The default action is to remove them")