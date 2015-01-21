/*
 * Copyright 2013-2014 Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Last modified on 15.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.sample.metadata

import org.apache.commons.io.FilenameUtils

/**
 * Some useful utils for metadata manipulation 
 */
public class MetadataUtil {
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
        defaultMetadataTable.createRow(sampleId, new ArrayList<String>())
    }

    /**
     * Gets a generic metadata table
     * @return metadata table which contains all statically-created metadata entries
     */
    public static MetadataTable getDefaultMetadataTable() {
        MetadataTable.GENERIC_METADATA_TABLE
    }
}
