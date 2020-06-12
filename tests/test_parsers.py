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

from pathlib import Path

from chunked_scatter.chunked_scatter import BedRegion, file_to_regions

datadir = Path(__file__).parent / Path("data")


def test_file_to_regions_bed():
    result = list(file_to_regions(datadir / "regions.bed"))
    assert result == [
        BedRegion("chr1", 100, 1000),
        BedRegion("chr1", 2000, 16000),
        BedRegion("chr2", 5000, 10000)
    ]


def test_file_to_regions_dict():
    result = list(file_to_regions(datadir / "ref.dict"))
    assert result == [
        BedRegion("chr1", 0, 3000000),
        BedRegion("chr2", 0, 500000)
    ]
