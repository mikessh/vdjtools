package com.antigenomics.vdjtools.graph;

import com.antigenomics.vdjtools.ClonotypeWrapper;
import com.antigenomics.vdjtools.ClonotypeWrapperContainer;
import com.antigenomics.vdjtools.join.ClonotypeKeyGen;
import com.antigenomics.vdjtools.overlap.OverlapType;
import com.antigenomics.vdjtools.sample.Clonotype;
import com.milaboratory.core.sequence.AminoAcidSequence;
import com.milaboratory.core.tree.NeighborhoodIterator;
import com.milaboratory.core.tree.SequenceTreeMap;

import java.util.*;
import java.util.function.Consumer;
import java.util.stream.StreamSupport;

public class DegreeStatisticsCalculator {
    private final ClonotypeGroupingFactory clonotypeGroupingFactory;
    private final GroupingSummary groupingSummary = new GroupingSummary();
    private final SequenceTreeMap<AminoAcidSequence, List<Clonotype>> stm = new SequenceTreeMap<>(AminoAcidSequence.ALPHABET);
    private final int substitutionThreshold, indelThreshold, totalMismatchThreshold;
    private final static ClonotypeKeyGen ckg = new ClonotypeKeyGen(OverlapType.Strict);

    public DegreeStatisticsCalculator(int substitutionThreshold,
                                      int indelThreshold,
                                      int totalMismatchThreshold,
                                      ClonotypeGroupingFactory clonotypeGroupingFactory) {
        this.substitutionThreshold = substitutionThreshold;
        this.indelThreshold = indelThreshold;
        this.totalMismatchThreshold = totalMismatchThreshold;
        this.clonotypeGroupingFactory = clonotypeGroupingFactory;
    }

    public <T extends ClonotypeWrapper> void inititalize(ClonotypeWrapperContainer<T> clonotypes) {
        final Map<AminoAcidSequence, List<Clonotype>> clonotypeMap = new HashMap<>();

        Spliterator<T> spliterator = Spliterators.spliterator(clonotypes.iterator(),
                clonotypes.getDiversity(),
                Spliterator.IMMUTABLE | Spliterator.SIZED);

        StreamSupport
                .stream(spliterator, true)
                .forEach(cw -> {
                    Clonotype clonotype = cw.getClonotype();
                    if (clonotype.isCoding()) {
                        clonotypeMap.computeIfAbsent(clonotype.getCdr3aaBinary(), k -> new ArrayList<>()).add(clonotype);
                        groupingSummary.update(clonotypeGroupingFactory.getGroup(clonotype)); // thread safe
                    }
                });

        for (Map.Entry<AminoAcidSequence, List<Clonotype>> entry : clonotypeMap.entrySet()) {
            stm.put(entry.getKey(), entry.getValue());
        }
    }

    public DegreeStatistics compute(Clonotype clonotype) {
        if (!clonotype.isCoding()) {
            return DegreeStatistics.UNDEF;
        }

        Set<Clonotype> clonotypeSet = new HashSet<>();
        ClonotypeGroup clonotypeGroup = clonotypeGroupingFactory.getGroup(clonotype);

        NeighborhoodIterator<AminoAcidSequence, List<Clonotype>> ni = stm.getNeighborhoodIterator(clonotype.getCdr3aaBinary(),
                substitutionThreshold, indelThreshold, indelThreshold, totalMismatchThreshold);
        List<Clonotype> matchList;

        while ((matchList = ni.next()) != null) {
            for (Clonotype match : matchList) {
                if (!ckg.generateKey(clonotype).equals(ckg.generateKey(match)) && // no match to itself
                        clonotypeGroup.equals(clonotypeGroupingFactory.getGroup(match))) { // same group, e.g. VJ length
                    clonotypeSet.add(match);
                }
            }
        }

        return new DegreeStatistics(clonotypeSet.size(),
                groupingSummary.getCount(clonotypeGroup));
    }
}
