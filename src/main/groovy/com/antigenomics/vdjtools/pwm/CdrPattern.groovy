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
 * Last modified on 23.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.pwm

import com.antigenomics.vdjtools.Clonotype

import java.util.regex.Pattern


class CdrPattern {
    private final int length
    private final String pattern
    private final Pattern _pattern
    private final double freq
    private final int count

    CdrPattern(String pattern, double freq, int count) {
        this.length = pattern.length()
        this.pattern = pattern
        this._pattern = Pattern.compile(pattern.replaceAll("X", "."))
        this.freq = freq
        this.count = count
    }

    int getLength() {
        length
    }

    String getPattern() {
        pattern
    }

    double getFreq() {
        return freq
    }

    int getCount() {
        return count
    }

    boolean matches(Clonotype clonotype) {
        _pattern.matcher(clonotype.cdr3aa).matches()
    }
}
