package com.antigenomics.vdjtools.db

/**
 * Created by mikesh on 10/23/14.
 */
class BrowserResult implements Iterable<CdrDatabaseMatch> {
    public final int sampleDiversity, databaseDiversity
    public final double sampleFreq
    private final List<CdrDatabaseMatch> matches

    BrowserResult(int sampleDiversity, int databaseDiversity, double sampleFreq, List<CdrDatabaseMatch> matches) {
        this.sampleDiversity = sampleDiversity
        this.databaseDiversity = databaseDiversity
        this.sampleFreq = sampleFreq
        this.matches = matches
    }

    int size() {
        matches.size()
    }

    @Override
    Iterator<CdrDatabaseMatch> iterator() {
        matches.iterator()
    }

    public static
    final String HEADER = "match_size\tsample_diversity_in_matches\tdb_diversity_in_matches\tsample_freq_in_matches"

    @Override
    String toString() {
        [size(), sampleDiversity, databaseDiversity, sampleFreq].join("\t")
    }
}
