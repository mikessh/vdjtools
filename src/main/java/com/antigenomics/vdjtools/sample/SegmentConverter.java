package com.antigenomics.vdjtools.sample;

import com.antigenomics.vdjtools.misc.Segment;
import com.antigenomics.vdjtools.misc.SegmentFactory;

import java.util.HashMap;
import java.util.Map;

public class SegmentConverter implements ClonotypeConverter {
    private final Map<String, Segment> vSegmentMap = new HashMap<>(),
            jSegmentMap = new HashMap<>();

    public SegmentConverter(Map<String, String> vSegmentMap, Map<String, String> jSegmentMap) {
        for (Map.Entry<String, String> conv : vSegmentMap.entrySet()) {
            this.vSegmentMap.put(conv.getKey(),
                    SegmentFactory.INSTANCE.create(conv.getValue()));
        }

        for (Map.Entry<String, String> conv : jSegmentMap.entrySet()) {
            this.jSegmentMap.put(conv.getKey(),
                    SegmentFactory.INSTANCE.create(conv.getValue()));
        }
    }

    @Override
    public Clonotype convert(Clonotype clonotype) {
        return clonotype.withSegments(
                vSegmentMap.getOrDefault(clonotype.getV(), clonotype.getVBinary()),
                jSegmentMap.getOrDefault(clonotype.getJ(), clonotype.getJBinary()));
    }
}
