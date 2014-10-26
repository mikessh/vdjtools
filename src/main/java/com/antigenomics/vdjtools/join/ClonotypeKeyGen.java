package com.antigenomics.vdjtools.join;

import com.antigenomics.vdjtools.Clonotype;
import com.antigenomics.vdjtools.intersection.IntersectionType;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

/**
 * Created by mikesh on 10/26/14.
 */
class ClonotypeKeyGen {
    private final IntersectionType intersectionType;

    public ClonotypeKeyGen(IntersectionType intersectionType) {
        this.intersectionType = intersectionType;
    }

    ClonotypeKey generateKey(Clonotype clonotype) {
        switch (intersectionType) {
            case Nucleotide:
                return new ClonotypeKey(clonotype) {
                    @Override
                    public boolean equals(Clonotype other) {
                        return clonotype.getCdr3nt().equals(other.getCdr3nt());
                    }

                    @Override
                    public int hashCode() {
                        return clonotype.getCdr3nt().hashCode();
                    }
                };

            case NucleotideV:
                return new ClonotypeKey(clonotype) {
                    @Override
                    public boolean equals(Clonotype other) {
                        return clonotype.getCdr3nt().equals(other.getCdr3nt()) &&
                                clonotype.getV().equals(other.getV());
                    }

                    @Override
                    public int hashCode() {
                        return clonotype.getCdr3nt().hashCode() * 31 + clonotype.getV().hashCode();
                    }
                };

            case NucleotideVJ:
                return new ClonotypeKey(clonotype) {
                    @Override
                    public boolean equals(Clonotype other) {
                        return clonotype.getCdr3nt().equals(other.getCdr3nt()) &&
                                clonotype.getV().equals(other.getV()) &&
                                clonotype.getJ().equals(other.getJ());
                    }

                    @Override
                    public int hashCode() {
                        return 31 * (clonotype.getCdr3nt().hashCode() * 31 + clonotype.getV().hashCode()) +
                                clonotype.getJ().hashCode();
                    }
                };

            case AminoAcid:
                return new ClonotypeKey(clonotype) {
                    @Override
                    public boolean equals(Clonotype other) {
                        return clonotype.getCdr3aa().equals(other.getCdr3aa());
                    }

                    @Override
                    public int hashCode() {
                        return clonotype.getCdr3aa().hashCode();
                    }
                };

            case AminoAcidV:
                return new ClonotypeKey(clonotype) {
                    @Override
                    public boolean equals(Clonotype other) {
                        return clonotype.getCdr3aa().equals(other.getCdr3aa()) &&
                                clonotype.getV().equals(other.getV());
                    }

                    @Override
                    public int hashCode() {
                        return clonotype.getCdr3aa().hashCode() * 31 + clonotype.getV().hashCode();
                    }
                };

            case AminoAcidVJ:
                return new ClonotypeKey(clonotype) {
                    @Override
                    public boolean equals(Clonotype other) {
                        return clonotype.getCdr3aa().equals(other.getCdr3aa()) &&
                                clonotype.getV().equals(other.getV()) &&
                                clonotype.getJ().equals(other.getJ());
                    }

                    @Override
                    public int hashCode() {
                        return 31 * (clonotype.getCdr3aa().hashCode() * 31 + clonotype.getV().hashCode()) +
                                clonotype.getJ().hashCode();
                    }
                };

            case AminoAcidNonNucleotide:
                return new ClonotypeKey(clonotype) {
                    @Override
                    public boolean equals(Clonotype other) {
                        return clonotype.getCdr3aa().equals(other.getCdr3aa()) &&
                                !clonotype.getCdr3nt().equals(other.getCdr3nt());
                    }

                    @Override
                    public int hashCode() {
                        return clonotype.getCdr3aa().hashCode();
                    }
                };

            case Strict:
                return new ClonotypeKey(clonotype) {
                    @Override
                    public boolean equals(Clonotype other) {
                        return clonotype.getKey().equals(other.getKey());
                    }

                    @Override
                    public int hashCode() {
                        return clonotype.getKey().hashCode();
                    }
                };

            default:
                throw new NotImplementedException();
        }
    }
}
