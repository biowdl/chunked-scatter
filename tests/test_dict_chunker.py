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

from chunked_scatter.chunked_scatter import dict_chunker


datadir = Path(__file__).parent / Path("data")


def test_dict_chunker():
    chunks = dict_chunker((datadir / Path("ref.dict")).open("r"), 1e6, 150)
    expected_output = [["chr1", 0, 1e6], ["chr1", 999850, 2e6],
                       ["chr1", 1999850, 3e6], ["chr2", 0, 5e5]]
    assert list(chunks) == expected_output


def test_dict_chunker_no_overlap():
    chunks = dict_chunker((datadir / Path("ref.dict")).open("r"), 1e6, 0)
    expected_output = [["chr1", 0, 1e6], ["chr1", 1e6, 2e6],
                       ["chr1", 2e6, 3e6], ["chr2", 0, 5e5]]
    assert list(chunks) == expected_output


def test_dict_chunker_big_value():
    chunks = dict_chunker((datadir / Path("ref.dict")).open("r"), 1e12, 0)
    expected_output = [["chr1", 0, 3e6], ["chr2", 0, 5e5]]
    assert list(chunks) == expected_output