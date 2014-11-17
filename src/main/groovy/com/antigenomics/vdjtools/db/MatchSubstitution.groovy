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
 * Last modified on 23.10.2014 by mikesh
 */

package com.antigenomics.vdjtools.db

import com.antigenomics.vdjtools.Clonotype

class MatchSubstitution {
    public final Clonotype parent
    public final int pos
    public final char to

    MatchSubstitution(Clonotype parent, int pos, char to) {
        this.parent = parent
        this.pos = pos
        this.to = to
    }

    public char getFrom() {
        parent.cdr3aa.charAt(pos)
    }

    @Override
    String toString() {
        "$pos:$from>$to"
    }
}
