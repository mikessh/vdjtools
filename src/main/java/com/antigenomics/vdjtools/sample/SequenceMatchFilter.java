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

package com.antigenomics.vdjtools.sample;

import com.antigenomics.vdjtools.misc.CommonUtil;
import org.apache.commons.lang3.StringUtils;

import java.util.regex.Pattern;

/**
 * A filter based on CDR3 pattern matching.
 */
public class SequenceMatchFilter extends ClonotypeFilter {
    private final static Pattern aaAccepted = Pattern.compile("[FLSYCWPHQRIMTNKVADEGX\\*_\\[\\]]+"),
            ntAccepted = Pattern.compile("[ATGCN\\[\\]]+");

    private final String patternString;
    private final Pattern pattern;
    private final boolean aminoAcid;
    private final int maxMismatches;

    public SequenceMatchFilter(String patternString, boolean aminoAcid, int maxMismatches) {
        super(false);
        patternString = patternString.toUpperCase();

        if (!(aminoAcid ? aaAccepted : ntAccepted).matcher(patternString).matches()) {
            throw new IllegalArgumentException("Bad sequence pattern: " + patternString);
        }

        if (StringUtils.countMatches(patternString, '[') !=
                StringUtils.countMatches(patternString, ']')) {
            throw new IllegalArgumentException("Bad sequence pattern: " + patternString);
        }

        patternString = patternString.replaceAll("[XN]", CommonUtil.PLACEHOLDER);

        this.pattern = Pattern.compile(patternString);
        this.patternString = patternString;
        this.aminoAcid = aminoAcid;
        this.maxMismatches = maxMismatches;
    }

    @Override
    protected boolean checkPass(Clonotype clonotype) {
        String query = aminoAcid ? clonotype.getCdr3aa() : clonotype.getCdr3nt();

        return pattern.matcher(query).find();
    }

    public String getPatternString() {
        return patternString;
    }

    public boolean isAminoAcid() {
        return aminoAcid;
    }

    public int getMaxMismatches() {
        return maxMismatches;
    }
}
