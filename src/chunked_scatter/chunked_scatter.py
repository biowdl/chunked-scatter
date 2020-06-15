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
from typing import Generator, Iterable, List, NamedTuple, Optional, Union


class BedRegion(NamedTuple):
    """A class that contains a region described as in the BED file format."""
    contig: str
    start: int
    end: int

    def __str__(self):
        return f"{self.contig}\t{self.start}\t{self.end}"


def dict_to_regions(in_file: Union[str, os.PathLike]
                    ) -> Generator[BedRegion, None, None]:
    """
    Converts a Picard SequenceDictionary file to a BedRegion Generator.
    :param in_file: The sequence dictionary
    :return: A generator of BedRegions
    """
    with open(in_file, "rt") as in_file_h:
        for line in in_file_h:
            fields = line.strip().split()
            if fields[0] != "@SQ":
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
    """
    Converts a BED file to a generator of BED regions
    :param in_file: The BED file
    :return: A BedRegion Generator
    """
    with open(in_file, "rt") as in_file_h:
        for line in in_file_h:
            fields = line.strip().split()
            # Skip browser and track fields and other invalid lines.
            if fields[0] in ["browser", "track"] or len(fields) < 3:
                continue
            # Take the first 3 columns of each line to create a new BedRegion
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
    Converts each region into chunks if the chunk_size is smaller than the
    region size.
    :param regions: The regions which to chunk.
    :param chunk_size: The size of the chunks.
    :param overlap: The size of the overlap between chunks.
    :return: The new chunked regions.
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
                    split_contigs: bool = False,
                    ) -> Generator[List[BedRegion], None, None]:
    """
    Scatter regions in chunks with an overlap. It returns Lists of regions
    where each list of regions has the regions describe at least the amount
    of base pairs in minimum bas pairs. Except the last list.
    :param regions: The regions which to chunk.
    :param chunk_size: The size of each chunk.
    :param overlap: How much overlap there should be between chunks.
    :param minimum_base_pairs: What the minimum amount of base pairs should be
    that the regions encompass per List.
    :param split_contigs: Whether contigs (chr1, for example) are allowed to be
    split across multiple lists.
    :return: Lists of BedRegions, which can be converted into BED files.
    """
    current_scatter_size = 0
    current_contig = None
    chunk_list: List[BedRegion] = []
    for chunk in region_chunker(regions, chunk_size, overlap):
        # If the next chunk is on a different contig
        if chunk.contig != current_contig or split_contigs:
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


def region_lists_to_scatter_files(region_lists: Iterable[List[BedRegion]],
                                  prefix: str) -> None:
    """
    Convert lists of BedRegions to '{prefix}{number}.bed' files. The number
    starts at 0 and is increased with 1 for each file.
    :param region_lists: The region lists to be converted into BED files.
    :param prefix: The filename prefix for the BedFiles
    :return: None
    """
    parent_dir = Path(prefix).parent
    if not parent_dir.exists():
        parent_dir.mkdir(parents=True)
    for scatter_number, region_list in enumerate(region_lists):
        out_file = f"{prefix}{scatter_number}.bed"
        with open(out_file, "wt") as out_file_h:
            for bed_region in region_list:
                out_file_h.write(str(bed_region) + os.linesep)


def main():
    args = parse_args()
    scattered_chunks = chunked_scatter(file_to_regions(args.input),
                                       args.chunk_size, args.overlap,
                                       args.minimum_bp_per_file,
                                       args.split_contigs)
    region_lists_to_scatter_files(scattered_chunks, args.prefix)


def common_parser() -> argparse.ArgumentParser:
    """Commmon arguments for chunked-scatter and scatter-regions."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--prefix", type=str, required=True,
                        help="The prefix of the ouput files. Output will be "
                        "named like: <PREFIX><N>.bed, in which N is an "
                        "incrementing number. This option is mandatory.")
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="The input file, either a bed file or a sequence "
                        "dict. Which format is used is detected by the "
                        "extension: '.bed' or '.dict'. This option is "
                        "mandatory.")
    parser.add_argument("--split-contigs", action="store_true",
                        help="If set, contigs are allowed to be split.")
    return parser


def parse_args():
    """Argument parser for the chunked-scatter program."""
    parser = common_parser()
    parser.description = (
        "Given a sequence dict or a bed file, scatter over the "
        "defined contigs/regions. Each contig/region will be split into "
        "multiple overlapping regions, which will be written to a new bed "
        "file. Each contig will be placed in a new file, unless the length of "
        "the contigs/regions doesn't exceed a given number.")

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


if __name__ == "__main__":  # pragma: no cover
    main()
