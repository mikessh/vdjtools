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

/**
 * An input stream factory. An example implementation would be creating a new file connection each time.
 * The main purpose of this class is to allow loading the same sample several times and to provide a mean
 * for implementing re-load for custom input streams.
 */
public interface InputStreamFactory {
    /**
     * Creates the corresponding input stream. Depending on stream type re-opens file connection, etc
     * @return a new input stream
     */
    public InputStream create()

    /**
     * Get the id associated with a given input stream
     * @return id string, a derivative of file name, etc
     */
    public String getId()
}