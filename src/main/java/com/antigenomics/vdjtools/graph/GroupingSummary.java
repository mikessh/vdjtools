package com.antigenomics.vdjtools.graph;


import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

public class GroupingSummary {
    private final Map<ClonotypeGroup, AtomicLong> counters = new ConcurrentHashMap<>();

    public void update(ClonotypeGroup clonotypeGroup) {
        counters.computeIfAbsent(clonotypeGroup, k -> new AtomicLong()).incrementAndGet();
    }

    public long getCount(ClonotypeGroup clonotypeGroup) {
        return counters.getOrDefault(clonotypeGroup, new AtomicLong()).get();
    }
}
