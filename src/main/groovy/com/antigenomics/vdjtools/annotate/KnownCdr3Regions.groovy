/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */

package com.antigenomics.vdjtools.annotate

import com.antigenomics.vdjtools.annotate.partitioning.Cdr3Center
import com.antigenomics.vdjtools.annotate.partitioning.Cdr3Region
import com.antigenomics.vdjtools.annotate.partitioning.DGermline
import com.antigenomics.vdjtools.annotate.partitioning.DJJunction
import com.antigenomics.vdjtools.annotate.partitioning.FullCdr3
import com.antigenomics.vdjtools.annotate.partitioning.JGermline
import com.antigenomics.vdjtools.annotate.partitioning.VDJunction
import com.antigenomics.vdjtools.annotate.partitioning.VGermline
import com.antigenomics.vdjtools.annotate.partitioning.VJJunction

class KnownCdr3Regions {
    private final Map<String, Cdr3Region> regionsByName

    static final KnownCdr3Regions INSTANCE = new KnownCdr3Regions()

    private KnownCdr3Regions() {
        this.regionsByName = [new VGermline(), new DGermline(), new JGermline(),
                              new VDJunction(), new DJJunction(),
                              new VJJunction(), new FullCdr3(),
                              new Cdr3Center()].collectEntries {
            [(it.name.toLowerCase()): it]
        }
    }

    Cdr3Region getByName(String name) {
        name = name.toLowerCase()
        if (!regionsByName.containsKey(name))
            throw new IllegalArgumentException("Bad CDR3 region name '$name', allowed values: ${allowedNames}")
        regionsByName[name]
    }

    List<String> getAllowedNames() {
        regionsByName.keySet().collect()
    }

    List<Cdr3Region> getAll() {
        regionsByName.values().collect()
    }
}
