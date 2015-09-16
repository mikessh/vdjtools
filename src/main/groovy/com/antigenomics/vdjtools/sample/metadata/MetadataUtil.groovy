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

package com.antigenomics.vdjtools.sample.metadata

import org.apache.commons.io.FilenameUtils

/**
 * Some useful utils for metadata manipulation 
 */
public class MetadataUtil {
    private static final Map<String, Integer> sampleHash = new HashMap<>()

    /**
     * Converts a file name to sample id 
     * @param fileName file name to convert
     * @return sample id, a shortcut for file name without any path and extension
     */
    public static String fileName2id(String fileName) {
        FilenameUtils.getBaseName(
                fileName.endsWith(".gz") ?
                        FilenameUtils.getBaseName(fileName) :
                        fileName)
    }

    /**
     * Creates sample metadata object and assigns it to a generic metadata table 
     * @param sampleId short unique identifier of a sample
     * @return sample metadata object assigned to a generic metadata table
     */
    public static SampleMetadata createSampleMetadata(String sampleId) {
        def idCount = (sampleHash[sampleId] ?: 0) + 1
        sampleHash.put(sampleId, idCount)
        defaultMetadataTable.createRow((idCount > 0 ? "$idCount." : "") + sampleId, new ArrayList<String>())
    }

    /**
     * Gets a generic metadata table
     * @return metadata table which contains all statically-created metadata entries
     */
    public static MetadataTable getDefaultMetadataTable() {
        MetadataTable.GENERIC_METADATA_TABLE
    }
}
