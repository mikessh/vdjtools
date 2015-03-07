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

import com.antigenomics.vdjtools.sample.Sample

/**
 * A sample provider
 */
public interface SampleConnection {
    /**
     * Gets the underlying sample. Note that in some implementations this could be very time consuming, 
     * as it will involve storing large amount of objects in memory and reading a file.
     * @return a sample object filled with clonotypes.
     */
    public Sample getSample()

    /**
     * Runs through sample file/stream without loading the whole sample into memory and collects statistics.
     * @return a blank sample object holding general statistics (clonotype count, etc).
     */
    public Sample haveAGlance()
}