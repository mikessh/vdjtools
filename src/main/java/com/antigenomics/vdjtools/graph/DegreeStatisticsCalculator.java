package com.antigenomics.vdjtools.graph;

import com.antigenomics.vdjtools.ClonotypeWrapper;
import com.antigenomics.vdjtools.ClonotypeWrapperContainer;
import com.antigenomics.vdjtools.sample.Clonotype;
import com.milaboratory.core.sequence.AminoAcidSequence;
import com.milaboratory.core.tree.NeighborhoodIterator;
import com.milaboratory.core.tree.SequenceTreeMap;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.StreamSupport;

public class DegreeStatisticsCalculator {
    private final ClonotypeGroupingFactory primaryClonotypeGroupingFactory,
            secondaryClonotypeGroupingFactory;
    private final GroupingSummary primaryGroupingSummary = new GroupingSummary(),
            secondaryGroupingSummary = new GroupingSummary();
    private final SequenceTreeMap<AminoAcidSequence, Queue<Clonotype>> stm = new SequenceTreeMap<>(AminoAcidSequence.ALPHABET);
    private final int substitutionThreshold, indelThreshold, totalMismatchThreshold;

    public DegreeStatisticsCalculator(int substitutionThreshold,
                                      int indelThreshold,
                                      int totalMismatchThreshold,
                                      ClonotypeGroupingFactory primaryClonotypeGroupingFactory) {
        this(substitutionThreshold, indelThreshold, totalMismatchThreshold,
                primaryClonotypeGroupingFactory, DummyClonotypeGroupingFactory.INSTANCE);
    }

    public DegreeStatisticsCalculator(int substitutionThreshold,
                                      int indelThreshold,
                                      int totalMismatchThreshold,
                                      ClonotypeGroupingFactory primaryClonotypeGroupingFactory,
                                      ClonotypeGroupingFactory secondaryClonotypeGroupingFactory) {
        this.substitutionThreshold = substitutionThreshold;
        this.indelThreshold = indelThreshold;
        this.totalMismatchThreshold = totalMismatchThreshold;
        this.primaryClonotypeGroupingFactory = primaryClonotypeGroupingFactory;
        this.secondaryClonotypeGroupingFactory = secondaryClonotypeGroupingFactory;
    }

    public <T extends ClonotypeWrapper> void inititalize(ClonotypeWrapperContainer<T> clonotypes) {
        final Map<AminoAcidSequence, Queue<Clonotype>> clonotypeMap = new ConcurrentHashMap<>();

        Spliterator<T> spliterator = Spliterators.spliterator(clonotypes.iterator(),
                clonotypes.getDiversity(),
                Spliterator.IMMUTABLE | Spliterator.SIZED);

        StreamSupport
                .stream(spliterator, true)
                .forEach(cw -> {
                    Clonotype clonotype = cw.getClonotype();
                    if (clonotype.isCoding()) {
                        clonotypeMap.computeIfAbsent(clonotype.getCdr3aaBinary(),
                                k -> new ConcurrentLinkedQueue<>()).add(clonotype);
                        primaryGroupingSummary.update(primaryClonotypeGroupingFactory.getGroup(clonotype)); // thread safe
                        secondaryGroupingSummary.update(secondaryClonotypeGroupingFactory.getGroup(clonotype));
                    }
                });

        for (Map.Entry<AminoAcidSequence, Queue<Clonotype>> entry : clonotypeMap.entrySet()) {
            stm.put(entry.getKey(), entry.getValue());
        }
    }

    public DegreeStatistics compute(Clonotype clonotype) {
        if (!clonotype.isCoding()) {
            return DegreeStatistics.UNDEF;
        }

        Set<Clonotype> clonotypeSet = new HashSet<>();
        ClonotypeGroup primaryClonotypeGroup = primaryClonotypeGroupingFactory.getGroup(clonotype),
                secondaryClonotypeGroup = secondaryClonotypeGroupingFactory.getGroup(clonotype);

        NeighborhoodIterator<AminoAcidSequence, Queue<Clonotype>> ni = stm.getNeighborhoodIterator(clonotype.getCdr3aaBinary(),
                substitutionThreshold, indelThreshold, indelThreshold, totalMismatchThreshold);
        Queue<Clonotype> matchList;

        while ((matchList = ni.next()) != null) {
            // explicit check for indel sum threshold as insertions and deletions are counted separately
            if (ni.getCurrentMutations().countOfIndels() <= indelThreshold) {
                for (Clonotype match : matchList) {
                    if (primaryClonotypeGroup.equals(primaryClonotypeGroupingFactory.getGroup(match))) { // same group, e.g. VJ length
                        clonotypeSet.add(match);
                    }
                }
            }
        }

        return new DegreeStatistics(clonotypeSet.size(),
                primaryGroupingSummary.getCount(primaryClonotypeGroup),
                secondaryGroupingSummary.getCount(secondaryClonotypeGroup));
    }
}
