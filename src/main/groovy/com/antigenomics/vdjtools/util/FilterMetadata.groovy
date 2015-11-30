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

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataEntryExpressionFilter

def cli = new CliBuilder(usage: "FilterMetadata [options] metadata.txt output_prefix")
cli.h("display help message")
cli.f(longOpt: "filter", argName: "string", args: 1, required: true,
        "Filter expression, metadata column names should be marked with ${MetadataEntryExpressionFilter.FILTER_MARK}, " +
                "e.g. \"__chain__=~/TR[AB]/\" or \"__chain__=='TRA'||__chain__=='TRB'\"")

def opt = cli.parse(args)

if (opt == null) {
    System.exit(2)
}

if (opt.h || opt.arguments().size() != 2) {
    cli.usage()
    System.exit(2)
}

def metadataFileName = opt.arguments()[0], filter = (String) opt.f, outputPrefix = opt.arguments()[1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

// Lazy load sample list, need to get absolute paths
println "[${new Date()} $scriptName] Checking sample(s)"
def sampleCollection = new SampleCollection((String) metadataFileName, Software.VDJtools, false)

println "[${new Date()} $scriptName] Filtering metadata by $filter"
def filteredMetadataTable = sampleCollection.metadataTable.select(new MetadataEntryExpressionFilter(filter))
filteredMetadataTable.storeWithOutput(outputPrefix, sampleCollection)
println "[${new Date()} $scriptName] Finished"