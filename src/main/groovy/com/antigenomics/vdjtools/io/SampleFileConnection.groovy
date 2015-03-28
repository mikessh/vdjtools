/*
 * Copyright 2013-2015 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 */

package com.antigenomics.vdjtools.io

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.metadata.SampleMetadata

/**
 * A wrapper for plain-text clonotype table stored in a file.
 * This is a semi-internal class to provide lazy-loading support for SampleCollection
 */
public class SampleFileConnection extends SampleStreamConnection {
    /**
     * Loads a sample from the specified file.
     * @param fileName path to file containing the sample.
     * @return sample object filled with clonotypes.
     */
    public static Sample load(String fileName) {
        load(fileName, Software.VDJtools)
    }

    /**
     * Loads a sample from the specified file.
     * @param fileName path to file containing the sample.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     * @return sample object filled with clonotypes.
     */
    public static Sample load(String fileName, Software software) {
        load(new FileInputStreamFactory(fileName), software)
    }

    /**
     * Loads a sample from the specified file.
     * @param fileName path to file containing the sample.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     * @param sampleMetadata a metadata object that will be associated with a given sample.
     * @return sample object filled with clonotypes.
     */
    public static Sample load(String fileName, Software software, SampleMetadata sampleMetadata) {
        load(new FileInputStreamFactory(fileName), software, sampleMetadata)
    }

    /**
     * Creates a sample connection, an object that could be used to access (load to memory, store, etc) a sample stored as plain text.
     * Will load the sample upon initialization and store it into memory. Generic sample metadata will be associated with the underlying sample.
     * @param fileName path to file containing the sample.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     */
    public SampleFileConnection(String fileName, Software software) {
        super(new FileInputStreamFactory(fileName), software)
    }

    /**
     * Creates a sample connection, an object that could be used to access (load to memory, store, etc) a sample stored as plain text.
     * Will load the sample upon initialization and store it into memory.
     * @param fileName path to file containing the sample.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     * @param sampleMetadata a metadata object that will be associated with a given sample.
     */
    public SampleFileConnection(String fileName, Software software, SampleMetadata sampleMetadata) {
        this(fileName, software, sampleMetadata, false, true)
    }

    /**
     * Creates a sample connection, an object that could be used to access (load to memory, store, etc) a sample stored as plain text.
     * @param fileName path to file containing the sample.
     * @param software type of software used to create clonotype table. Specifies how the plain-text input will be parsed.
     * @param sampleMetadata a metadata object that will be associated with a given sample.
     * @param lazy will load the sample only when {@code getSample ( )} is called.
     * @param store sample will be stored into memory after loading.
     */
    public SampleFileConnection(String fileName, Software software, SampleMetadata sampleMetadata,
                                boolean lazy, boolean store) {
        super(new FileInputStreamFactory(fileName), software, sampleMetadata, lazy, store)
    }
}
