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

package com.antigenomics.vdjtools.util.db

import com.antigenomics.vdjdb.core.db.CdrEntry
import com.antigenomics.vdjdb.core.query.CdrSearchResult
import com.antigenomics.vdjtools.sample.Clonotype
import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence

/**
 * A class representing a match between CDR3 database entry and a clonotype   
 */
public class CdrMatch {
    private final CdrEntry subject
    private final boolean vMatch, jMatch
    private final Clonotype query
    private final Alignment<AminoAcidSequence> alignment

    public static List<CdrMatch> collectMatches(Clonotype query,
                                                CdrSearchResult cdrSearchResult,
                                                boolean requireVMatch, boolean requireJMatch) {
        def result = new ArrayList<CdrMatch>()
        cdrSearchResult.cdrEntrySet.each {
            def wrapper = new CdrMatch(it, query, cdrSearchResult.alignment)
            if ((!requireVMatch || wrapper.vMatch) && (!requireJMatch || wrapper.jMatch))
                result.add(wrapper)
        }
        result
    }

    public CdrMatch(CdrEntry subject, Clonotype query,
                    Alignment<AminoAcidSequence> alignment) {
        this.subject = subject
        this.query = query
        this.vMatch = query.v == subject.v //todo: re-implement V stuff
        this.jMatch = query.j == subject.j
        this.alignment = alignment
    }

    /**
     * Gets the CDR3 database entry 
     * @return
     */
    public CdrEntry getSubject() {
        subject
    }

    /**
     * Gets the clonotype 
     * @return
     */
    public Clonotype getQuery() {
        query
    }

    /**
     * Returns @code{true} if the Variable segment of the clonotype is the same as in CDR3 database entry, otherwise @{code false}* @return
     */
    public boolean isvMatch() {
        vMatch
    }

    /**
     * Returns @code{true} if the Joining segment of the clonotype is the same as in CDR3 database entry, otherwise @{code false}* @return
     */
    public boolean isjMatch() {
        jMatch
    }

    /**
     * Gets the alignment between clonotype's CDR3 amino acid sequence and the matching CDR3 from database 
     * @return
     */
    public Alignment<AminoAcidSequence> getAlignment() {
        alignment
    }

    /**
     * Header string, used for tabular output
     */
    public static final String HEADER = "score\tquery_cdr3aa\tquery_v\tquery_j\t" +
            "subject_cdr3aa\tsubject_v\tsubject_j\t" +
            "v_match\tj_match\tmismatches"

    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        [alignment.score,
         query.cdr3aa, query.v, query.j,
         subject.cdr3aa, subject.v, subject.j,
         vMatch, jMatch, alignment.absoluteMutations.toString()
        ].flatten().join("\t")
    }
}
