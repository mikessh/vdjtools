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

import com.antigenomics.vdjtools.sample.metadata.MetadataUtil
import com.antigenomics.vdjtools.util.CommonUtil

/**
 * A file input stream factory. This factory creates a new file connection each time.
 */
public class FileInputStreamFactory implements InputStreamFactory {
    private final String fileName

    /**
     * Creates a new instance of file input stream factory associated with a given file name
     * @param fileName path to underlying file
     */
    public FileInputStreamFactory(String fileName) {
        this.fileName = fileName
    }

    /**
     * @inheritDoc
     */
    @Override
    public InputStream create() {
        CommonUtil.getFileStream(fileName)
    }

    /**
     * @inheritDoc
     */
    @Override
    public String getId() {
        MetadataUtil.fileName2id(fileName)
    }
}
