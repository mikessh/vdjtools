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

package com.antigenomics.vdjtools.io

import com.antigenomics.vdjtools.misc.Software
import org.junit.Test

import static Software.*
import static com.antigenomics.vdjtools.TestUtil.getResource
import static com.antigenomics.vdjtools.io.SampleStreamConnection.load

class SampleLoadTest {
    private static void loadTest(Software software, int count, int diversity) {
        loadTest(software, '', count, diversity);
    }

    private static void loadTest(Software software, String suffix, int count, int diversity) {
        def resStream = getResource("samples/${software.toString().toLowerCase()}${suffix}.txt.gz")
        def sample = load(resStream, software)

        assert sample.count == count
        assert sample.diversity == diversity
        if (software.perReadOutput) {
            // check if frequency is recalculated
            assert Math.abs(sample.freqAsInInput - 1.0) < 1e-5
        }
    }

    @Test
    public void mitcrTest() {
        loadTest(MiTcr, 10000, 6493)
    }

    @Test
    public void migecTest() {
        loadTest(MiGec, 6147, 2420)
    }

    @Test
    public void migmapTest() {
        loadTest(MigMap, 720, 703)
    }

    @Test
    public void simpleTest() {
        loadTest(VDJtools, 11878, 2257)
    }

    @Test
    public void immunoseqTest() {
        loadTest(ImmunoSeq, 2345835, 13052)
    }

    @Test
    public void imgtTest() {
        loadTest(ImgtHighVQuest, 9681, 7199)
    }

    @Test
    public void mixcrTest() {
        loadTest(MiXcr, 96132, 262)
    }

    @Test
    public void mixcrFullLengthTest() {
        loadTest(MiXcr, ".fl", 14156, 33)
    }

    @Test
    public void mixcr17LengthTest() {
        loadTest(MiXcr, ".171", 900, 873)
    }

    @Test
    public void imseqTest() {
        loadTest(ImSeq, 647, 86)
    }

    @Test
    public void emptyTest() {
        def resStream = getResource("samples/empty.txt")
        def sample = load(resStream)

        assert sample.count == 0
        assert sample.diversity == 0
    }
}
