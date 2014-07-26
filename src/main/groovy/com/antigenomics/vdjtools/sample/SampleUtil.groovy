/**
 Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

package com.antigenomics.vdjtools.sample

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.Software
import org.apache.commons.io.FilenameUtils
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation

class SampleUtil {
    static List<Clonotype> loadClonotypes(String fileName, Software software) {
        def clonotypes = new ArrayList()
        def inputFile = new File(fileName)
        inputFile.withReader { reader ->
            for (int i = 0; i < software.headerLineCount; i++)
                reader.readLine()

            def line
            while ((line = reader.readLine()) != null) {
                if (!software.comment || !line.startsWith(software.comment))
                    clonotypes.add(Clonotype.parseClonotype(line, software))
            }
        }
        clonotypes
    }

    static Sample loadSample(String fileName, Software software) {
        new Sample(new SampleMetadata(FilenameUtils.getBaseName(fileName)),
                loadClonotypes(fileName, software))
    }

    static double correlation(List<Clonotype> clonotypes1, List<Clonotype> clonotypes2) {
        if (clonotypes1.size() != clonotypes2.size())
            throw new IllegalArgumentException("Clonotype set sizes should match")

        if (clonotypes1.size() < 2)
            return Double.NaN

        def freq1 = clonotypes1.collect { Math.log10(it.freq) } as double[],
            freq2 = clonotypes2.collect { Math.log10(it.freq) } as double[]

        new PearsonsCorrelation().correlation(freq1, freq2)
    }
}
