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

package com.antigenomics.vdjtools.annotate.partitioning;

import com.antigenomics.vdjtools.sample.Clonotype;
import com.milaboratory.core.Range;

public class Cdr3Center extends Cdr3Region {
    private final int span;

    public static final Cdr3Center CDR3_CENTER_5 = new Cdr3Center(2),  // extracts -2,-1,0,+1,+2 amino acids
            CDR3_CENTER_3 = new Cdr3Center(1);  // extracts -1,0,+1 amino acids

    public Cdr3Center(int span) {
        this.span = span;
    }

    @Override
    protected Range getRange(Clonotype clonotype) {
        int cdr3CenterAA = clonotype.getCdr3Length() / 6;
        return new Range(Math.max(0, 3 * (cdr3CenterAA - span)),
                Math.min(clonotype.getCdr3Length(), 3 * (cdr3CenterAA + span + 1)));
    }

    public int getSpan() {
        return span;
    }

    @Override
    public String getName() {
        return "CDR3-center-" + (span * 2 + 1);
    }
}
