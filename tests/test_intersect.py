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

from chunked_scatter.chunked_scatter import BedRegion, intersect_bed_regions

def test_intersect():
    regions_a = [
        BedRegion("chr1", 0, 14000),
        BedRegion("chr2", 0, 15000),
        BedRegion("chr3", 5000, 16000)
    ]
    regions_b = [
        BedRegion("chr0", 0, 5000),
        BedRegion("chr1", 3000, 5000),
        BedRegion("chr1", 7000, 13000),
        BedRegion("chr1", 13500, 16000),
        BedRegion("chr1", 17000, 21000),
        BedRegion("chr2", 10000, 19000),
        BedRegion("chr3", 0, 10000)
    ]
    result = list(intersect_bed_regions(regions_a, regions_b))
    assert result == [
        BedRegion("chr1", 3000, 5000),
        BedRegion("chr1", 7000, 13000),
        BedRegion("chr1", 13500, 14000),
        BedRegion("chr2", 10000, 15000),
        BedRegion("chr3", 5000, 10000)
    ]