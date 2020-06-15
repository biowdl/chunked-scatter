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
from typing import Generator, Iterable, List

from .chunked_scatter import chunked_scatter, common_parser, \
    region_lists_to_scatter_files
from .parsers import BedRegion, file_to_regions

DEFAULT_SCATTER_SIZE = 10**9


def merge_regions(regions: Iterable[BedRegion]
                  ) -> Generator[BedRegion, None, None]:
    """
    Merge regions that overlap or are exactly adjacent
    :param regions: An iterable of possibly overlapping regions
    :return: a generator of merged regions
    """
    merged_region = None
    for region in regions:
        if merged_region is None:
            merged_region = region
        else:
            if (merged_region.contig == region.contig and
                    merged_region.end >= region.start and
                    region.end >= merged_region.start):
                start = min(merged_region.start, region.start)
                end = max(merged_region.end, region.end)
                merged_region = BedRegion(merged_region.contig, start, end)
            else:
                yield merged_region
                merged_region = region
    if merged_region:
        yield merged_region


def scatter_regions(regions: Iterable[BedRegion],
                    scattersize: int,
                    contigs_can_be_split: bool = False,
                    ) -> Generator[List[BedRegion], None, None]:
    """
    Interface to chunked_scatter with sane defaults that make it function
    similar to biopet-scatterregions. Also any overlapping bed records are
    merged together.
    :param regions: The regions over which to scatter.
    :param scattersize: What the size of the scatter should be.
    :param contigs_can_be_split: Whether contigs (chr1, for example) are
    allowed to be split across multiple lists.
    :return: Yields lists of BedRegions which can be converted into bed files.
    """
    region_lists = chunked_scatter(regions,
                                   chunk_size=scattersize,
                                   minimum_base_pairs=scattersize,
                                   overlap=0,
                                   contigs_can_be_split=contigs_can_be_split)
    for region_list in region_lists:
        yield list(merge_regions(region_list))


def argument_parser() -> argparse.ArgumentParser:
    """Argument parser for the scatter-regions program."""
    parser = common_parser()
    parser.description = (
        "Given a sequence dict or a bed file, scatter over the defined "
        "contigs/regions. Creates a bed file where the contigs add up "
        "approximately to the given scatter size.")
    parser.add_argument("-s", "--scatter-size", type=int,
                        default=DEFAULT_SCATTER_SIZE,
                        help=f"How large the regions over which to scatter "
                             f"should be. Default: {DEFAULT_SCATTER_SIZE}.")
    return parser


def main():
    args = argument_parser().parse_args()
    scattered_chunks = scatter_regions(file_to_regions(args.input),
                                       args.scatter_size,
                                       contigs_can_be_split=args.split_contigs)
    region_lists_to_scatter_files(scattered_chunks, args.prefix)
