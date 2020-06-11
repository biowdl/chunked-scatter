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
import sys

import pytest

from chunked_scatter.chunked_scatter import parse_args, main, input_is_bed


datadir = Path(__file__).parent / Path("data")


def test_bed_input(tmpdir):
    sys.argv = ["script", "-p", "{}/test_result_".format(tmpdir), "-i",
                str(Path(datadir, "regions.bed")), "-c", "5000", "-m", "1"]
    main()
    file1 = tmpdir / Path("test_result_0.bed")
    file2 = tmpdir / Path("test_result_1.bed")
    expected_file1 = ("chr1\t100\t1000\n"
                      "chr1\t2000\t7000\n"
                      "chr1\t6850\t12000\n"
                      "chr1\t11850\t16000\n")
    expected_file2 = "chr2\t5000\t10000\n"
    print(file1, file1.exists())
    print(file2, file2.exists())
    assert file1.exists()
    assert file2.exists()
    assert file1.read() == expected_file1
    assert file2.read() == expected_file2


def test_dict_input(tmpdir):
    sys.argv = ["script", "-p", "{}/test_result_".format(tmpdir), "-i",
                str(Path(datadir, "ref.dict"))]
    main()
    file1 = tmpdir / Path("test_result_0.bed")
    expected_file1 = ("chr1\t0\t1000000\n"
                      "chr1\t999850\t2000000\n"
                      "chr1\t1999850\t3000000\n"
                      "chr2\t0\t500000\n")
    print(file1, file1.exists())
    assert file1.exists()
    assert file1.read() == expected_file1


def test_check_input_extension_wrong_ext(capsys):
    with pytest.raises(ValueError,
        match="Only files with .bed or .dict extensions are supported."):
        input_is_bed(Path("input"))


def test_check_input_extension_bed():
    assert input_is_bed(Path("input.bed"))


def test_check_input_extension_dict():
    assert not input_is_bed(Path("input.dict"))


def test_parse_args():
    sys.argv = ["script", "-p", "prefix", "-i", "input.bed"]
    args = parse_args()
    assert args.input == Path("input.bed")
    assert args.prefix == "prefix"
    assert args.chunk_size == 1e6
    assert args.minimum_bp_per_file == 45e6
    assert args.overlap == 150
