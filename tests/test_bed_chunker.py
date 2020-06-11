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

from chunked_scatter.chunked_scatter import bed_file_to_regions, region_chunker


datadir = Path(__file__).parent / Path("data")


def test_bed_chunker():
    chunks = region_chunker(bed_file_to_regions(Path(datadir, "regions.bed")
                                                ), 5000, 150)
    expected_output = [["chr1", 100, 1000], ["chr1", 2000, 7000],
                       ["chr1", 6850, 12000], ["chr1", 11850, 16000],
                       ["chr2", 5000, 10000]]
    assert [list(chunk) for chunk in chunks] == expected_output


def test_bed_chunker_no_overlap():
    chunks = region_chunker(bed_file_to_regions(Path(datadir, "regions.bed")
                                                ), 5000, 0)
    expected_output = [["chr1", 100, 1000], ["chr1", 2000, 7000],
                       ["chr1", 7000, 12000], ["chr1", 12000, 16000],
                       ["chr2", 5000, 10000]]
    assert [list(chunk) for chunk in chunks] == expected_output