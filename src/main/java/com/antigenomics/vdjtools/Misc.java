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
 * Last modified on 30.1.2015 by mikesh
 */

package com.antigenomics.vdjtools;

/**
 * Some miscellaneous constants and functions
 */
public class Misc {
    /**
     * Upper limit on current precision of RepSeq
     */
    public static final double JITTER = 1e-9, JITTER_LOG10 = Math.log10(JITTER);
}
