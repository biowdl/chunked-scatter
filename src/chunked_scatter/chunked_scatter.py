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
from typing import TextIO


def dict_chunker(in_file: TextIO, chunk_size: int, overlap: int):
    """
    Chunk the chromosomes/contigs from a dict.
    :param in_file: The input file.
    :param chunk_size: The size of the chunks.
    :param overlap: The size of the overlap between chunks.
    """
    for line in in_file.readlines():
        fields = line.split()
        if fields[0] == "@SQ":
            for field in fields:
                if field.startswith("LN"):
                    length = int(field[3:])
                elif field.startswith("SN"):
                    name = field[3:]
            # This will cause the last chunk to be between 0.5 and 1.5
            # times the chunk_size in length, this way we avoid the
            # possibility that the last chunk ends up being to small
            # (eg. 1+overlap bases).
            position = 0
            while position + chunk_size*1.5 < length:
                if position-overlap <= 0:
                    yield [name, 0, int(position+chunk_size)]
                else:
                    yield [name, int(position-overlap),
                           int(position+chunk_size)]
                position += chunk_size
            if position-overlap <= 0:
                yield [name, 0, length]
            else:
                yield [name, int(position-overlap), length]


def bed_chunker(in_file: TextIO, chunk_size: int, overlap: int):
    """
    Chunk the entries of the bed file.
    :param in_file: The input file.
    :param chunk_size: The size of the chunks.
    :param overlap: The size of the overlap between chunks.
    """
    for line in in_file.readlines():
        fields = line.strip().split("\t")
        if fields[0] not in ["browser", "track"] and len(line) >= 3:
            start = int(fields[1])
            end = int(fields[2])
            name = fields[0]
            position = start
            # This will cause the last chunk to be between 0.5 and 1.5
            # times the chunk_size in length, this way we avoid the
            # possibility that the last chunk ends up being to small
            # (eg. 1+overlap bases).
            while position + chunk_size*1.5 < end:
                if position-overlap <= start:
                    yield [name, start, int(position+chunk_size)]
                else:
                    yield [name, int(position-overlap),
                           int(position+chunk_size)]
                position += chunk_size
            if position-overlap <= start:
                yield [name, start, end]
            else:
                yield [name, int(position-overlap), end]


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
    with args.input.open("r") as in_file:
        if bed_input:
            chunks = bed_chunker(in_file, args.chunk_size, args.overlap)
        else:
            chunks = dict_chunker(in_file, args.chunk_size, args.overlap)

        for chunk in chunks:
            # If the next chunk is on a different contig
            if chunk[0] != current_contig:
                current_contig = chunk[0]
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
            current_contents += "{}\t{}\t{}\n".format(chunk[0], chunk[1],
                                                      chunk[2])
            current_scatter_size += (chunk[2]-chunk[1])
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
