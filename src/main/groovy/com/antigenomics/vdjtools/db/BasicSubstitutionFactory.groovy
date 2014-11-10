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
 * Last modified on 10.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.db

import com.antigenomics.vdjtools.Clonotype

public class BasicSubstitutionFactory implements SubstitutionFactory {
    @Override
    public MatchSubstitution create(Clonotype parent, int pos, char from, char to) {
        (to != from &&
                pos > 0 && pos < (parent.cdr3aa.length() - 1)) ? // conserved Cys & Phe
                new MatchSubstitution(parent, pos, to) : null
    }
}