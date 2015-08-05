/*
 * Copyright 2013-2015 Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.antigenomics.vdjtools;

import java.util.HashMap;
import java.util.Map;

public class SegmentFactory {
    public static final SegmentFactory INSTANCE = new SegmentFactory();

    protected final Map<String, Segment> segmentCache = new HashMap<>();

    public SegmentFactory() {
        segmentCache.put(Segment.MISSING.name, Segment.MISSING);
    }

    public Segment create(String name) {
        Segment segment = segmentCache.get(name);

        if (segment == null) {
            segmentCache.put(name, segment = new Segment(name));
        }

        return segment;
    }
}
