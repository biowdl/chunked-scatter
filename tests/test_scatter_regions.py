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

from chunked_scatter.chunked_scatter import BedRegion
from chunked_scatter.scatter_regions import merge_regions, scatter_regions

import pytest

from .test_chunkers import DICT_REGIONS


MERGE_REGION_TESTS = [
 ([BedRegion("chr1", 0, 10000), BedRegion("chr1", 10000, 20000)],
  [BedRegion("chr1", 0, 20000)]),
 ([BedRegion("chr1", 0, 9999), BedRegion("chr1", 10000, 20000)],
  [BedRegion("chr1", 0, 9999), BedRegion("chr1", 10000, 20000)]),
 ([BedRegion("chr1", 0, 10000), BedRegion("chr1", 8000, 20000)],
  [BedRegion("chr1", 0, 20000)])
]


@pytest.mark.parametrize(["regions", "result"], MERGE_REGION_TESTS)
def test_adjacent_regions(regions, result):
    merged = list(merge_regions(regions))
    assert merged == result


def test_scatter_regions():
    result = list(scatter_regions(DICT_REGIONS, 1_000_000))
    assert result == [
        [BedRegion("chr1", 0, 3_000_000)],
        [BedRegion("chr2", 0, 500_000)]
    ]


def test_scatter_regions_split_contigs():
    result = list(scatter_regions(DICT_REGIONS, 900_000, split_contigs=True))
    assert result == [
        [BedRegion("chr1", 0, 900_000)],
        [BedRegion("chr1", 900_000, 1_800_000)],
        [BedRegion("chr1", 1_800_000, 2_700_000)],
        [BedRegion("chr1", 2_700_000, 3_000_000),
         BedRegion("chr2", 0, 500_000)]
    ]
