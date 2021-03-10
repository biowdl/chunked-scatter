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

import sys
from pathlib import Path

from chunked_scatter.chunked_scatter import main, parse_args
from chunked_scatter.safe_scatter import main as safe_scatter_main
from chunked_scatter.scatter_regions import main as scatter_regions_main

DATA_DIR = Path(__file__).parent / Path("data")


def test_bed_input(tmpdir, capsys):
    sys.argv = ["script", "-p", "{}/test_result_".format(tmpdir),
                str(Path(DATA_DIR, "regions.bed")), "-c", "5000", "-m", "1",
                "-P"]
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
    assert str(tmpdir / Path("test_result_0.bed")) in capsys.readouterr().out


def test_dict_input(tmpdir):
    sys.argv = ["script", "-p", "{}/test_result_".format(tmpdir),
                str(Path(DATA_DIR, "ref.dict"))]
    main()
    file1 = tmpdir / Path("test_result_0.bed")
    expected_file1 = ("chr1\t0\t1000000\n"
                      "chr1\t999850\t2000000\n"
                      "chr1\t1999850\t3000000\n"
                      "chr2\t0\t500000\n")
    print(file1, file1.exists())
    assert file1.exists()
    assert file1.read() == expected_file1


def test_parse_args():
    sys.argv = ["script", "input.bed"]
    args = parse_args()
    assert args.input == "input.bed"
    assert args.prefix == "scatter-"
    assert args.chunk_size == 1e6
    assert args.minimum_bp_per_file == 45e6
    assert args.overlap == 150


def test_scatter_regions_main(tmpdir, capsys):
    sys.argv = ["scatter-regions", "-p",
                str(Path(str(tmpdir), "scatters", "scatter-")),
                "--split-contigs",
                "-s", "1100000",
                str(Path(DATA_DIR, "ref.dict")),
                "--print-paths"]
    scatter_regions_main()
    assert Path(str(tmpdir), "scatters", "scatter-0.bed").exists()
    assert Path(str(tmpdir), "scatters", "scatter-1.bed").exists()
    assert Path(str(tmpdir), "scatters", "scatter-2.bed").exists()
    assert Path(str(tmpdir), "scatters", "scatter-3.bed").exists()
    assert Path(str(tmpdir), "scatters", "scatter-2.bed").read_text() == (
        "chr1\t2200000\t3000000\n"
    )
    captured = capsys.readouterr()
    assert str(Path(str(tmpdir), "scatters", "scatter-0.bed")) in captured.out


def test_safe_scatter_main(tmpdir, capsys):
    sys.argv = ["safe-scatter", "-p",
                str(Path(str(tmpdir), "scatters", "scatter-")),
                "--scatter-count", "3",
                str(Path(DATA_DIR, "ref.dict")),
                "--print-paths"]
    safe_scatter_main()
    assert Path(str(tmpdir), "scatters", "scatter-0.bed").exists()
    assert Path(str(tmpdir), "scatters", "scatter-1.bed").exists()
    assert Path(str(tmpdir), "scatters", "scatter-2.bed").exists()
    assert not Path(str(tmpdir), "scatters", "scatter-3.bed").exists()
    assert Path(str(tmpdir), "scatters", "scatter-1.bed").read_text() == (
        "chr1\t1160000\t2320000\n"
    )
    captured = capsys.readouterr()
    assert str(Path(str(tmpdir), "scatters", "scatter-0.bed")) in captured.out
