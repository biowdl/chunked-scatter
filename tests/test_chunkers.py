# Copyright (c) 2019 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import pytest

from chunked_scatter.chunked_scatter import BedRegion, chunked_scatter, \
    region_chunker

BED_REGIONS = [BedRegion("chr1", 100, 1000), BedRegion("chr1", 2000, 16000),
               BedRegion("chr2", 5000, 10000)]

DICT_REGIONS = [BedRegion("chr1", 0, 3000000), BedRegion("chr2", 0, 500000)]

REGION_TESTS = [
    (BED_REGIONS, 5000, 150, [
        BedRegion("chr1", 100, 1000),
        BedRegion("chr1", 2000, 7000),
        BedRegion("chr1", 6850, 12000),
        BedRegion("chr1", 11850, 16000),
        BedRegion("chr2", 5000, 10000)
    ]),
    (BED_REGIONS, 5000, 0, [
        BedRegion("chr1", 100, 1000),
        BedRegion("chr1", 2000, 7000),
        BedRegion("chr1", 7000, 12000),
        BedRegion("chr1", 12000, 16000),
        BedRegion("chr2", 5000, 10000)
    ]),
    (DICT_REGIONS, 1_000_000, 150, [
        BedRegion("chr1", 0, 1_000_000),
        BedRegion("chr1", 999850, 2_000_000),
        BedRegion("chr1", 1999850, 3_000_000),
        BedRegion("chr2", 0, 500_000)
    ]),
    (DICT_REGIONS, 1_000_000, 0, [
        BedRegion("chr1", 0, 1_000_000),
        BedRegion("chr1", 1_000_000, 2_000_000),
        BedRegion("chr1", 2_000_000, 3_000_000),
        BedRegion("chr2", 0, 500_000)
    ]),
    (DICT_REGIONS, 10**12, 150, [
        BedRegion("chr1", 0, 3_000_000),
        BedRegion("chr2", 0, 500_000)
    ])
]


@pytest.mark.parametrize(["regions", "chunk_size", "overlap", "result"],
                         REGION_TESTS)
def test_region_chunker(regions, chunk_size, overlap, result):
    chunks = list(region_chunker(regions, chunk_size, overlap))
    assert chunks == result
