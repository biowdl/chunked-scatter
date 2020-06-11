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

import argparse
from pathlib import Path
from typing import NamedTuple, Generator, Optional, Iterable


class BedRegion(NamedTuple):
    contig: str
    start: int
    end: int

    def __str__(self):
        return f"{self.contig}\t{self.start}\t{self.end}"


def dict_to_regions(in_file: str) -> Generator[BedRegion, None, None]:
    with open(in_file, "rt") as in_file_h:
        for line in in_file_h:
            fields = line.strip().split()
            if not fields[0] == "@SQ":
                continue

            contig: Optional[str] = None
            length: Optional[int] = None
            for field in fields:
                if field.startswith("LN"):
                    length = int(field[3:])
                elif field.startswith("SN"):
                    contig = field[3:]
            if contig and length:
                yield BedRegion(contig, 0, length)


def bed_file_to_regions(in_file: str) -> Generator[BedRegion, None, None]:
    with open(in_file, "rt") as in_file_h:
        for line in in_file_h:
            # Take the first 3 columns of each line to create a new BedRegion
            fields = line.strip().split()
            yield BedRegion(fields[0], int(fields[1]), int(fields[2]))


def region_chunker(regions: Iterable[BedRegion], chunk_size: int, overlap: int
                   ) -> Generator[BedRegion, None, None]:
    """
    :param regions:
    :param chunk_size: The size of the chunks.
    :param overlap: The size of the overlap between chunks.
    :return:
    """
    for contig, start, end in regions:
        position = start
        # This will cause the last chunk to be between 0.5 and 1.5
        # times the chunk_size in length, this way we avoid the
        # possibility that the last chunk ends up being to small
        # (eg. 1+overlap bases).
        while position + chunk_size * 1.5 < end:
            if position - overlap <= start:
                yield BedRegion(contig, start, int(position + chunk_size))
            else:
                yield BedRegion(contig, int(position - overlap),
                                int(position + chunk_size))
            position += chunk_size
        if position - overlap <= start:
            yield BedRegion(contig, start, end)
        else:
            yield BedRegion(contig, int(position - overlap), end)


def input_is_bed(filename: Path):
    """
    Check whether or not the input file is a bed file.
    :param filename: The filename.
    """
    if filename.suffix == ".bed":
        return True
    elif filename.suffix == ".dict":
        return False
    else:
        raise ValueError(
            "Only files with .bed or .dict extensions are supported.")


def main():
    args = parse_args()
    bed_input = input_is_bed(args.input)
    current_scatter = 0
    current_scatter_size = 0
    current_contig = None
    current_contents = str()
    if bed_input:
        regions = bed_file_to_regions(args.input)
    else:
        regions = dict_to_regions(args.input)
    chunked_regions = region_chunker(regions, args.chunk_size, args.overlap)

    for chunk in chunked_regions:
        # If the next chunk is on a different contig
        if chunk.contig != current_contig:
            current_contig = chunk.contig
            # and the current bed file contains enough bases
            if current_scatter_size >= args.minimum_bp_per_file:
                # write the bed file
                with open("{}{}.bed".format(args.prefix, current_scatter),
                          "w") as out_file:
                    out_file.write(current_contents)
                # and start compiling the next bed file
                current_scatter += 1
                current_scatter_size = 0
                current_contents = str()
        # Add the chunk to the bed file
        current_contents += "{}\t{}\t{}\n".format(*chunk)
        current_scatter_size += (chunk.end-chunk.start)
    # Write last bed file
    with open("{}{}.bed".format(args.prefix, current_scatter),
              "w") as out_file:
        out_file.write(current_contents)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Given a sequence dict or a bed file, scatter over the "
        "defined contigs/regions. Each contig/region will be split into "
        "multiple overlapping regions, which will be written to a new bed "
        "file. Each contig will be placed in a new file, unless the length of "
        "the contigs/regions doesn't exceed a given number.")
    parser.add_argument("-p", "--prefix", type=str, required=True,
                        help="The prefix of the ouput files. Output will be "
                        "named like: <PREFIX><N>.bed, in which N is an "
                        "incrementing number. This option is mandatory.")
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="The input file, either a bed file or a sequence "
                        "dict. Which format is used is detected by the "
                        "extension: '.bed' or '.dict'. This option is "
                        "mandatory.")
    parser.add_argument("-c", "--chunk-size", type=int, default=1e6,
                        metavar="SIZE",
                        help="The size of the chunks. The first chunk in a "
                        "region or contig will be exactly length SIZE, "
                        "subsequent chunks will SIZE + OVERLAP and the final "
                        "chunk may be anywhere from 0.5 to 1.5 times SIZE "
                        "plus overlap. If a region (or contig) is smaller "
                        "than SIZE the original regions will be returned. "
                        "Defaults to 1e6")
    parser.add_argument("-m", "--minimum-bp-per-file", type=int, default=45e6,
                        help="The minimum number of bases represented within "
                        "a single output bed file. If an input contig or "
                        "region is smaller than this MINIMUM_BP_PER_FILE, "
                        "then the next contigs/regions will be placed in the "
                        "same file untill this minimum is met. Defaults to "
                        "45e6.")
    parser.add_argument("-o", "--overlap", type=int, default=150,
                        help="The number of bases which each chunk should "
                        "overlap with the preceding one. Defaults to 150.")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
