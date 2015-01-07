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
 *
 * Last modified on 7.1.2015 by mikesh
 */

package com.antigenomics.vdjtools.io

import com.antigenomics.vdjtools.sample.metadata.MetadataUtil
import com.antigenomics.vdjtools.util.CommonUtil

class FileInputStreamFactory implements InputStreamFactory {
    private final String fileName

    public FileInputStreamFactory(String fileName) {
        this.fileName = fileName
    }

    @Override
    InputStream create() {
        CommonUtil.getFileStream(fileName)
    }

    @Override
    String getId() {
        MetadataUtil.fileName2id(fileName)
    }
}
