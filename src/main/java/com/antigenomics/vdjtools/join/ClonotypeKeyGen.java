package com.antigenomics.vdjtools.join;

import com.antigenomics.vdjtools.Clonotype;
import com.antigenomics.vdjtools.intersection.IntersectionType;
import com.antigenomics.vdjtools.join.key.*;
import com.antigenomics.vdjtools.sample.Sample;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by mikesh on 10/26/14.
 */
public class ClonotypeKeyGen {
    private final IntersectionType intersectionType;

    public ClonotypeKeyGen(IntersectionType intersectionType) {
        this.intersectionType = intersectionType;
    }

    public Set<ClonotypeKey> generateKeySet(Sample sample) {
        Set<ClonotypeKey> keySet = new HashSet<>();
        for (Clonotype clonotype : sample) {
            keySet.add(generateKey(clonotype));
        }
        return keySet;
    }

    public ClonotypeKey generateKey(Clonotype clonotype) {
        switch (intersectionType) {
            case Nucleotide:
                return new NtKey(clonotype);

            case NucleotideV:
                return new NtVKey(clonotype);

            case NucleotideVJ:
                return new NtVJKey(clonotype);

            case AminoAcid:
                return new AaKey(clonotype);

            case AminoAcidV:
                return new AaVKey(clonotype);

            case AminoAcidVJ:
                return new AaVJKey(clonotype);

            case AminoAcidNonNucleotide:
                return new AaNotNtKey(clonotype);

            case Strict:
                return new StrictKey(clonotype);

            default:
                throw new NotImplementedException();
        }
    }
}
