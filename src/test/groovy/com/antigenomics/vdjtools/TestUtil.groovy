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


package com.antigenomics.vdjtools

import com.antigenomics.vdjtools.io.InputStreamFactory
import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.sample.SampleCollection

import java.util.zip.GZIPInputStream

import static com.antigenomics.vdjtools.io.SampleStreamConnection.load

class TestUtil {
    static final SampleCollection DEFAULT_SAMPLE_COLLECTION = loadSamples(),
            SINGLE_EMPTY_SAMPLE = SampleCollection.fromSampleList([load(getResource("samples/empty.txt"), Software.VDJtools)])

    private static SampleCollection loadSamples() {
        def samples = Software.values().collect {
            load(getResource("samples/${it.toString().toLowerCase()}.txt.gz"), it)
        }
        samples.add(load(getResource("samples/empty.txt"), Software.VDJtools))

        SampleCollection.fromSampleList(samples)
    }


    public static InputStreamFactory getResource(String resourceName) {
        [
                create: {
                    def is = TestUtil.class.classLoader.getResourceAsStream(resourceName)
                    resourceName.endsWith(".gz") ? new GZIPInputStream(is) : is
                },
                getId : { resourceName.split("/")[-1] }
        ] as InputStreamFactory
    }
}
