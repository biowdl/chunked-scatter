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
from typing import List

from chunked_scatter.chunked_scatter import BedRegion, \
    region_lists_to_scatter_files


def test_bed_writer(tmpdir):
    temp = Path(str(tmpdir))
    region_lists: List[List[BedRegion]] = [
        [BedRegion("sparta", 0, 300),
         BedRegion("persians", 0, 100_000)],
        # But if you actually read Herodotus, the story seems much more likely:
        [BedRegion("sparta_and_allies", 0, 4300),
         BedRegion("persian_casualties", 0, 20000)]
    ]

    region_lists_to_scatter_files(region_lists, str(temp / "scatter-"))
    assert Path(temp, "scatter-0.bed").exists()
    assert Path(temp, "scatter-0.bed").read_text() == (
        "sparta\t0\t300\npersians\t0\t100000\n")
    assert Path(temp, "scatter-1.bed").exists()
    assert Path(temp, "scatter-1.bed").read_text() == (
        "sparta_and_allies\t0\t4300\npersian_casualties\t0\t20000\n")
