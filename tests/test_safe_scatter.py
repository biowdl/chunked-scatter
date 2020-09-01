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

from chunked_scatter import safe_scatter
from chunked_scatter.chunked_scatter import BedRegion

import pytest

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
BIN_SIZE_TESTS = [
    ([BedRegion("chr1", 0, 100)], 10, 10),
    ([BedRegion("chr1", 0, 100)], 3, 33),
    ([BedRegion("chr1", 0, 100)], 7, 14),
    ([BedRegion("chr1", 0, 100)], 9, 11),
    ([BedRegion("chr1", 0, 50), BedRegion("chr1", 50, 100)], 9, 11)
]

# Region, min_scatter_size, scattered regions
SCATTER_REGIONS_TESTS = [
    ([BedRegion("chr1", 0, 100)], 20, [
        BedRegion("chr1", 0, 20),
        BedRegion("chr1", 20, 40),
        BedRegion("chr1", 40, 60),
        BedRegion("chr1", 60, 80),
        BedRegion("chr1", 80, 100)]),
    ([BedRegion("chr1", 0, 100)], 27, [
        BedRegion("chr1", 0, 27),
        BedRegion("chr1", 27, 54),
        BedRegion("chr1", 54, 100)]),
    ([BedRegion("chr1", 0, 5)], 1, [
       BedRegion("chr1", 0, 1),
       BedRegion("chr1", 1, 2),
       BedRegion("chr1", 2, 3),
       BedRegion("chr1", 3, 4),
       BedRegion("chr1", 4, 5)]),
    ([BedRegion("chr1", 0, 2)], 1.2, [
       BedRegion("chr1", 0, 1),
       BedRegion("chr1", 1, 2)]),
    ([BedRegion("chr1", 0, 10)], 20, [
       BedRegion("chr1", 0, 10)]),
    ([BedRegion("chr1", 0, 10), BedRegion("chr2", 0, 12)], 20, [
       BedRegion("chr1", 0, 10), BedRegion("chr2", 0, 12)])
]

SCATTER_REGIONS_INVALID = [
    ([BedRegion("chr1", 0, 5)], 0),
    ([BedRegion("chr1", 0, 5)], -1),
    ([BedRegion("chr1", 0, 5)], 0.2)
]

# region, scatter_count, scatter_size
SAFE_SCATTER_INVALID = [
    ([BedRegion("chr1", 0, 100)], 4, 100),
    ([BedRegion("chr1", 0, 100)], 4, 26)
]

# region, scatter_count, scatter_size, result
SAFE_SCATTER_TESTS = [
    ([BedRegion("chr1", 0, 100)], 1, 50, [
        [BedRegion("chr1", 0, 100)]
    ]),
    ([BedRegion("chr1", 0, 100)], 1, 100, [
        [BedRegion("chr1", 0, 100)]
    ]),
    ([BedRegion("chr1", 0, 100)], 4, 25, [
        [BedRegion("chr1", 0, 25)],
        [BedRegion("chr1", 25, 50)],
        [BedRegion("chr1", 50, 75)],
        [BedRegion("chr1", 75, 100)]
    ])
]


@pytest.mark.parametrize(["regions", "result"], MERGE_REGION_TESTS)
def test_adjacent_regions(regions, result):
    merged = list(safe_scatter.merge_regions(regions))
    assert merged == result


@pytest.mark.parametrize(["regions", "result"], SUM_REGION_TESTS)
def test_sum_regions(regions, result):
    safe_scatter.sum_regions(regions) == result


@pytest.mark.parametrize(["regions", "scatter_count", "result"],
                         BIN_SIZE_TESTS)
def test_determine_bin_size(regions, scatter_count, result):
    assert safe_scatter.determine_bin_size(regions, scatter_count) == result


@pytest.mark.parametrize(["regions", "min_scatter_size", "result"],
                         SCATTER_REGIONS_TESTS)
def test_scatter_regions(regions, min_scatter_size, result):
    scatter_result = list(safe_scatter.scatter_regions(regions,
                          min_scatter_size))
    assert scatter_result == result


@pytest.mark.parametrize(["regions", "min_scatter_size"],
                         SCATTER_REGIONS_INVALID)
def test_scatter_regions_sanity(regions, min_scatter_size):
    with pytest.raises(RuntimeError):
        next(safe_scatter.scatter_regions(regions, min_scatter_size))


@pytest.mark.parametrize(["regions", "scatter_count", "min_scatter_size"],
                         SAFE_SCATTER_INVALID)
def test_safe_scatter_sanity(regions, scatter_count, min_scatter_size):
    with pytest.raises(RuntimeError):
        next(safe_scatter.safe_scatter(regions, scatter_count,
             min_scatter_size))


@pytest.mark.parametrize(["regions", "scatter_count", "min_scatter_size",
                          "result"], SAFE_SCATTER_TESTS)
def test_safe_scatter(regions, scatter_count, min_scatter_size, result):
    scattered_regions = list(safe_scatter.safe_scatter(regions, scatter_count,
                             min_scatter_size))
    assert scattered_regions == result
