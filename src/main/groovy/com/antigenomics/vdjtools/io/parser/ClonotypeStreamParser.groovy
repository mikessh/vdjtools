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


package com.antigenomics.vdjtools.io.parser

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample

/**
 * Base class for providing parsing of various RepSeq software output.
 * The stream parser is not thread-safe. 
 */
public abstract class ClonotypeStreamParser implements Iterable<Clonotype> {
    protected final List<String> header = new ArrayList<>()
    protected final Software software
    protected final Iterator<String> innerIter
    protected final Sample sample
    private int skippedLineCount = 0, commentLineCount
    private boolean printedWarning = false
    private final boolean hasComment

    /**
     * Creates a new instance of clonotype parser. 
     * It is a clonotype factory wrapped around the specified input stream
     * @param innerIter object that iterates over file lines, i.e. rows in the clonotype table
     * @param software software used to create the clonotype table. Specifies the parser
     * @param sample a blank sample to fill up with clonotypes
     */
    protected ClonotypeStreamParser(Iterator<String> innerIter, Software software, Sample sample) {
        this.software = software
        this.innerIter = innerIter
        this.sample = sample
        this.hasComment = software.comment && software.comment.length() > 0
        this.commentLineCount = software.headerLineCount
    }

    /**
     * Gets an instance of clonotype parser.
     * It is a clonotype factory wrapped around the specified input stream
     * @param inputStream input stream that will be read
     * @param software software used to create the clonotype table. Specifies the parser
     * @param sample a blank sample to fill up with clonotypes
     * @return clonotype parser object
     */
    public static ClonotypeStreamParser create(InputStream inputStream, Software software, Sample sample) {
        ClonotypeStreamParser parser
        def reader = new BufferedReader(new InputStreamReader(inputStream))
        def innerIter = reader.iterator()

        switch (software) {
            case Software.MiTcr:
                parser = new MiTcrParser(innerIter, sample)
                break
            case Software.IgBlast:
                parser = new IgBlastParser(innerIter, sample)
                break
            case Software.Simple:
                parser = new SimpleParser(innerIter, sample)
                break
            case Software.MiGec:
                parser = new MiGecParser(innerIter, sample)
                break
            case Software.ImmunoSeq:
                parser = new ImmunoSeqParser(innerIter, sample)
                break
            case Software.ImgtHighVQuest:
                parser = new ImgtHighVQuestParser(innerIter, sample)
                break
            default:
                throw new UnsupportedOperationException("Don't know how to parse $software data")
        }

        (0..<software.headerLineCount).each { parser.header.add(innerIter.next()) }

        return parser
    }

    /**
     * Parses a string into clonotype in a {@code Software}-dependent manner 
     * @param clonotypeString string to parse
     * @return string to parse
     */
    protected abstract Clonotype innerParse(String clonotypeString)

    /**
     * Parses a string into clonotype in a {@code Software}-dependent manner and performs some consistency checks.
     * Skips input strings that result in incomplete / bad clonotypes.
     * @param clonotypeString string to parse
     * @return a clonotype instance or {@code null} if input string was skipped
     */
    public Clonotype parse(String clonotypeString) {
        if (hasComment && clonotypeString.startsWith(software.comment)) {
            commentLineCount++
            return null
        }

        def clonotype = innerParse(clonotypeString)

        def badFieldMap = clonotype ? ["no_cdr3nt" : missingEntry(clonotype.cdr3nt),
                                       "no_cdr3aa" : missingEntry(clonotype.cdr3aa),
                                       "no_v"      : missingEntry(clonotype.v),
                                       "no_j"      : missingEntry(clonotype.j),
                                       "zero_count": clonotype.count == 0,
                                       "zero_freq" : !software.perReadOutput && clonotype.freqAsInInput == 0] :
                ["bad_line": true]

        if (badFieldMap.any { it.value }) {
            if (!printedWarning) {
                printedWarning = true
                println "[WARNING] Some of the essential fields are bad/missing " +
                        "for the following clonotype string (displaying first 5 warnings)"
            }
            if (skippedLineCount++ < 5) {
                println badFieldMap.findAll { it.value }.collect { it.key }.join(",") + ":"
                println "$clonotypeString"
            }
            return null
        }

        return clonotype
    }

    /**
     * Gets the number of skipped input lines
     * @return number of input lines that resulted in bad/incomplete clonotypes
     */
    public int getSkippedLineCount() {
        skippedLineCount
    }

    /**
     * Gets the number of comment input lines 
     * @return number of input lines that have a comment character at their beginning
     */
    public int getCommentLineCount() {
        commentLineCount
    }

    /**
     * As for now, just reports summary statistics to {@code stdout}
     */
    public void finish() {
        println "[${new Date()} ClonotypeStreamParser] Finished parsing. " +
                "$commentLineCount header and $skippedLineCount bad line(s) were skipped."
    }

    /**
     * INTERNAL checks if a given entry string is missing
     * @param entry
     * @return
     */
    private static boolean missingEntry(String entry) {
        !entry || entry.length() == 0 || entry == "."
    }

    /**
     * Implementation that simplifies usage syntax {@code parser.each{ Clonotype c -> ...}}
     * @return a clonotype iterator wrapped around the inner plain text table iterator
     */
    @Override
    public Iterator<Clonotype> iterator() {
        [hasNext: {
            innerIter.hasNext()
        }, next : {
            parse(innerIter.next())
        }] as Iterator
    }
}
