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
from chunked_scatter.safe_scatter import *

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

SUM_REGION_TESTS = [
    ([BedRegion("chr1", 0, 100)], 100),
    ([BedRegion("chr1", 0, 100), BedRegion("chr2", 0, 100)], 200),
    ([BedRegion("chr1", 0, 0)], 0)
]

# Region, scatter_count, determined scatter size
SCATTER_SIZE_TESTS = [
    ([BedRegion("chr1", 0, 100)], 10, 10)
]
@pytest.mark.parametrize(["regions", "result"], MERGE_REGION_TESTS)
def test_adjacent_regions(regions, result):
    merged = list(merge_regions(regions))
    assert merged == result


@pytest.mark.parametrize(["regions", "result"], SUM_REGION_TESTS)
def test_sum_regions(regions, result):
    sum_regions(regions) == result


@pytest.mark.parametrize(["regions", "scatter_count", "result"],
        SCATTER_SIZE_TESTS)
def test_determine_scatter_size(regions, scatter_count, result):
    assert determine_scatter_size(regions, scatter_count) == result
