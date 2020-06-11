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
import os
from pathlib import Path
from typing import NamedTuple, Generator, Optional, Iterable, List, Union


class BedRegion(NamedTuple):
    contig: str
    start: int
    end: int

    def __str__(self):
        return f"{self.contig}\t{self.start}\t{self.end}"


def intersect_bed_regions(regions_a: Iterable[BedRegion],
                          regions_b: Iterable[BedRegion]
                          ) -> Generator[BedRegion, None, None]:
    # We iterate over b_regions multiple times.
    b_regions = list(regions_b)

    for region_a in regions_a:
        for region_b in b_regions:
            if not region_a.contig == region_b.contig:
                continue
            if region_a.start < region_b.end and region_b.start < region_a.end:
                start = max(region_a.start, region_b.start)
                end = min(region_a.end, region_b.end)
                yield BedRegion(region_a.contig, start, end)


def dict_to_regions(in_file: Union[str, os.PathLike]
                    ) -> Generator[BedRegion, None, None]:
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


def bed_file_to_regions(in_file: Union[str, os.PathLike]
                        ) -> Generator[BedRegion, None, None]:
    with open(in_file, "rt") as in_file_h:
        for line in in_file_h:
            # Take the first 3 columns of each line to create a new BedRegion
            fields = line.strip().split()
            yield BedRegion(fields[0], int(fields[1]), int(fields[2]))


def file_to_regions(in_file: Path):
    if in_file.suffix == ".bed":
        return bed_file_to_regions(in_file)
    elif in_file.suffix == ".dict":
        return dict_to_regions(in_file)
    else:
        raise NotImplementedError(
            "Only files with .bed or .dict extensions are supported.")


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


def chunked_scatter(regions: Iterable[BedRegion],
                    chunk_size: int,
                    overlap: int,
                    minimum_base_pairs: int,
                    intersect_regions: Optional[Iterable[BedRegion]] = None
                    ) -> Generator[List[BedRegion], None, None]:
    if intersect_regions is not None:
        regions = intersect_bed_regions(regions, intersect_regions)
    current_scatter_size = 0
    current_contig = None
    chunk_list = []
    for chunk in region_chunker(regions, chunk_size, overlap):
        # If the next chunk is on a different contig
        if chunk.contig != current_contig:
            current_contig = chunk.contig
            # and the current bed file contains enough bases
            if current_scatter_size >= minimum_base_pairs:
                yield chunk_list
                chunk_list = []
                current_scatter_size = 0
        # Add the chunk to the bed file
        chunk_list.append(chunk)
        current_scatter_size += (chunk.end-chunk.start)
    # If there are leftovers yield them.
    if chunk_list:
        yield chunk_list


def main():
    args = parse_args()
    regions = file_to_regions(args.input)
    scattered_chunks = chunked_scatter(regions, args.chunk_size, args.overlap,
                            args.minimum_bp_per_file)
    for current_scatter, chunk_list in enumerate(scattered_chunks):
        out_file =  f"{args.prefix}{current_scatter}.bed"
        with open(out_file, "wt") as out_file_h:
            out_file_h.writelines(str(bed_regions) + "\n"
                                  for bed_regions in chunk_list)


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
