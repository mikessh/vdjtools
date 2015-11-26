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

package com.antigenomics.vdjtools.util

import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.BlankMetadataEntryFilter

def cli = new CliBuilder(usage: "SplitMetadata [options] metadata.txt output_prefix")
cli.h("display help message")
cli.c(longOpt: "columns", argName: "string1,string2,...", args: 1, required: true,
        "Column name(s) to split metadata by.")

def opt = cli.parse(args)

if (opt.h || opt.arguments().size() != 2) {
    cli.usage()
    System.exit(1)
}

def metadataFileName = opt.arguments()[0], columnIds = ((String) opt.c).split(","), outputPrefix = opt.arguments()[1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

// Lazy load sample list, need to get absolute paths
println "[${new Date()} $scriptName] Checking sample(s)"
def sampleCollection = new SampleCollection((String) metadataFileName)

println "[${new Date()} $scriptName] Splitting metadata by $columnIds"

def sampleIdByMetadataValue = new HashMap<String, List<String>>()

sampleCollection.each { sample ->
    def key = columnIds.collect { sample.sampleMetadata[it].value }.join(".")
    def sampleList = sampleIdByMetadataValue[key]
    if (sampleList == null) {
        sampleIdByMetadataValue.put(key, sampleList = new ArrayList<String>())
    }
    sampleList.add(sample.sampleMetadata.sampleId)
}

sampleIdByMetadataValue.each {
    def filteredMetadataTable = sampleCollection.metadataTable.select(BlankMetadataEntryFilter.INSTANCE,
            new HashSet<String>(it.value))
    filteredMetadataTable.storeWithOutput(outputPrefix, sampleCollection, it.key)
}

println "[${new Date()} $scriptName] Finished"
