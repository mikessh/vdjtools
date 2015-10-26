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

package com.antigenomics.vdjtools.sample;

import com.antigenomics.vdjtools.ClonotypeContainer;
import com.antigenomics.vdjtools.Software;
import com.antigenomics.vdjtools.io.parser.ClonotypeStreamParser;
import com.antigenomics.vdjtools.sample.metadata.SampleMetadata;

import java.io.InputStream;
import java.util.*;

public class Sample implements ClonotypeContainer {
    private final List<Clonotype> clonotypes = new ArrayList<>();
    private final SampleMetadata sampleMetadata;
    private double frequency = 0;
    private long count = 0;
    private int diversity = 0;

    private Sample(SampleMetadata sampleMetadata) {
        this.sampleMetadata = sampleMetadata;
    }

    public Sample(Sample toCopy, Map<Clonotype, Integer> samplerMap) {
        this.sampleMetadata = toCopy.sampleMetadata;

        for (Clonotype clonotype : toCopy.clonotypes) {
            Integer newCount = samplerMap.get(clonotype);

            if (newCount != null && newCount > 0)
                this.addClonotype(new Clonotype(clonotype, this, newCount));
        }

        Collections.sort(clonotypes);
    }

    public Sample(Sample toCopy, ClonotypeFilter filter, int top) {
        this.sampleMetadata = toCopy.sampleMetadata;

        for (Clonotype clonotype : toCopy.clonotypes) {
            if (top > -1 && this.getDiversity() == top)
                break;

            if (filter.pass(clonotype))
                this.addClonotype(new Clonotype(clonotype, this));
        }
    }

    public Sample(Sample toClone, ClonotypeFilter filter) {
        this(toClone, filter, -1);
    }

    public Sample(Sample toClone) {
        this(toClone, BlankClonotypeFilter.INSTANCE, -1);
    }

    public static Sample fromInputStream(InputStream inputStream,
                                         SampleMetadata sampleMetadata,
                                         Software software,
                                         int top, boolean store, boolean collapse) {
        Sample sample = new Sample(sampleMetadata);

        ClonotypeStreamParser clonotypeStreamParser = ClonotypeStreamParser.create(inputStream, software, sample);

        boolean sorted = !collapse;
        int prevCount = Integer.MAX_VALUE;

        Map<Clonotype, Clonotype> existingClonotypes = new HashMap<>();

        for (Clonotype clonotype : clonotypeStreamParser) {
            if (top > -1 && sample.getDiversity() == top)
                break;

            if (clonotype != null) {
                int count = clonotype.getCount();

                if (sorted && count > prevCount) {
                    sorted = false;
                    prevCount = count;
                }

                Clonotype existing = null;

                if (collapse) {
                    existing = existingClonotypes.get(clonotype);

                    if (existing != null) {
                        existing.append(clonotype);
                    } else {
                        existingClonotypes.put(clonotype, clonotype);
                    }
                }

                sample.addClonotype(clonotype, store, existing);
            }
        }

        clonotypeStreamParser.finish(); // report progress

        // on-demand sorting
        if (!sorted)
            Collections.sort(sample.clonotypes);

        // Re-calculate frequencies for per read storing software
        if (software.isPerReadOutput()) {
            sample.frequency = 0;
            for (Clonotype clonotype : sample) {
                sample.frequency += clonotype.recalculateFrequency();
            }
        }

        return sample;
    }

    public static Sample fromInputStream(InputStream inputStream,
                                         SampleMetadata sampleMetadata,
                                         Software software) {
        return fromInputStream(inputStream, sampleMetadata, software, -1, true, false);
    }

    private void addClonotype(Clonotype clonotype) {
        addClonotype(clonotype, true, null);
    }

    private void addClonotype(Clonotype clonotype, boolean store, Clonotype existingClonotype) {
        count += clonotype.getCount();
        frequency += clonotype.getFreqAsInInput();

        if (existingClonotype == null) {
            diversity++;
            if (store)
                clonotypes.add(clonotype);
        }
    }

    public SampleMetadata getSampleMetadata() {
        return sampleMetadata;
    }

    /**
     * For 1.8 stream
     * @return
     */
    public List<Clonotype> getClonotypes() {
        return Collections.unmodifiableList(clonotypes);
    }

    @Override
    public double getFreq() {
        return frequency;
    }

    @Override
    public long getCount() {
        return count;
    }

    @Override
    public int getDiversity() {
        return diversity;
    }

    @Override
    public Clonotype getAt(int index) {
        if (index < 0 || index >= clonotypes.size())
            throw new IndexOutOfBoundsException();
        return clonotypes.get(index);
    }

    @Override
    public boolean isSorted() {
        return true;
    }

    @Override
    public Iterator<Clonotype> iterator() {
        return clonotypes.iterator();
    }
}
