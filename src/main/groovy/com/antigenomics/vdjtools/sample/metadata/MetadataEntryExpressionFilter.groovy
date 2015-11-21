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

package com.antigenomics.vdjtools.sample.metadata

class MetadataEntryExpressionFilter implements MetadataEntryFilter {
    static final String FILTER_MARK = "__"
    final String expression
    final List<String> columnIds = new ArrayList<>()

    MetadataEntryExpressionFilter(String expression) {
        expression.split(FILTER_MARK).each { token ->
            def pattern = "$FILTER_MARK$token$FILTER_MARK"
            if (expression.contains(pattern)) {
                expression = expression.replaceAll(pattern, "x[\"$token\"]")
                columnIds.add(token)
            }
        }

        if (expression.contains(FILTER_MARK)) {
            throw new RuntimeException("Failed to parse filter, '$FILTER_MARK' symbols left.")
        }

        this.expression = expression
    }

    @Override
    boolean passes(Map<String, Object> entryValueMap) {
        if (!entryValueMap.keySet().containsAll(columnIds)) {
            throw new RuntimeException("Cannot process metadata $entryValueMap. " +
                    "Some of the $columnIds column IDs required for evaluating filter are missing.")
        }
        
        Eval.x(entryValueMap, expression)
    }
}
