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

from .chunked_scatter import common_parser, region_lists_to_scatter_files
from .parsers import BedRegion, file_to_regions


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


def sum_regions(regions: List[BedRegion]):
    """ Calculate the total length of all regions """
    return sum(len(region) for region in regions)


def determine_bin_size(regions: List[BedRegion],
                       scatter_count: int):
    """
    Determine the target scatter size, based on the total size of the regions
    and the number of target regions (scatter_count).
    """
    total_size = sum_regions(regions)
    return int(total_size/scatter_count)


def scatter_regions(regions: List[BedRegion],
                    min_scatter_size):
    """
    Scatter the regions into chunks. All chunks will be of size
    min_scatter_size, except (possibly) the last region.
    The last region will be >= min_scatter_size < min_scatter_size*2
    """
    # Make sure we don't get a floating minimum scatter size
    min_scatter_size = int(min_scatter_size)

    if min_scatter_size < 1:
        raise RuntimeError("min_scatter_size must be a positive integer")
    for region in regions:
        # How much of the region is still remaining
        remaining = len(region)

        # If the current region is so small we cannot get 2 min_scatter_size
        # from it, just yield the the entire region and continue
        if remaining < 2*min_scatter_size:
            yield region
            continue

        # We keep looping as long as there are at least 2 min_scatter_size left
        # of the region
        while remaining >= 2*min_scatter_size:
            contig, start, end = region
            new_region = BedRegion(contig, start, start+min_scatter_size)
            yield new_region
            region = BedRegion(contig, start+min_scatter_size, end)
            remaining = len(region)
        yield BedRegion(contig, start+min_scatter_size, end)


def safe_scatter(regions: List[BedRegion],
                 scatter_count: int,
                 min_scatter_size: int = 10000,
                 split_contigs: bool = True,
                 ) -> Generator[List[BedRegion], None, None]:
    """
    Scatter the regions equally over the specified scatter_count.

    The terminology here is that we will generate scatter_count "bins", and we
    will place each of the scattered regions (len >= min_scatter_size) in one
    of the bins.

    This is done in a manner that guarantees the following:
      - No regions will be split smaller than min_scatter_size. There can be
        regions smaller than min_scatter_size, but only when the original
        the original region was smaller than min_scatter_size.
      - All scatters will be within min_scatter_size of each other.

    Any overlapping bed records are merged together.
    :param regions: The regions over which to scatter.
    :param scatter_count: The number of bins to create.
    :param min_scatter_size: The minimum size of a scattered region.
    :param split_contigs: Whether contigs (chr1, for example) are
    allowed to be split across multiple lists.
    :return: Yields lists of BedRegions which can be converted into bed files.
    """
    # What is the target size for the bins?
    target_bin_size = determine_bin_size(regions, scatter_count)

    if target_bin_size < min_scatter_size:
        msg = (f"--min-scatter-size is not compatible with the provided "
               f"regions and number of bins ({min_scatter_size} > "
               f"{target_bin_size})")
        raise RuntimeError(msg)

    # Where we store all bins
    current_bin = None

    for region in scatter_regions(regions, min_scatter_size):
        # If this is the first ever region we parse, initialse the bin
        if current_bin is None:
            current_bin = [region]
            current_bin_size = len(region)
            continue
        # If adding this region would put us over the target bin size,
        # yield the bin and start a new one
        if current_bin_size + len(region) > target_bin_size:
            # Here we merge the chunks back together if they are adjacent
            yield list(merge_regions(current_bin))
            current_bin = [region]
            current_bin_size = len(region)
        # If this region does not put us over the target bin size, add it
        else:
            current_bin.append(region)
            current_bin_size += len(region)
    # If we are done, yield the last bin
    yield list(merge_regions(current_bin))


def argument_parser() -> argparse.ArgumentParser:
    """Argument parser for the scatter-regions program."""
    parser = common_parser()
    parser.description = (
        "Given a sequence dict, fasta index or a bed file, scatter over the "
        "defined contigs/regions. Creates a bed file where the contigs add up "
        "to the average scatter size to within min_scatter_size.")
    parser.add_argument("-c", "--scatter-count", type=int,
                        default=50,
                        help=f"The number of chunks to scatter the regions in."
                             f"All chunks will be within --min-scatter-size "
                             f"of each other.")
    parser.add_argument("-m", "--min-scatter-size", type=int,
                        default=10000,
                        help=f"The minimum size of a scatter. This tool will "
                             f"never generate regions smaller than this "
                             f"value, unless the original regions are"
                             f"smaller.")
    return parser


def main():
    args = argument_parser().parse_args()
    # We need all regions instead of an iterator
    regions = list(file_to_regions(args.input))
    scattered_chunks = list(safe_scatter(regions, args.scatter_count,
                                         args.min_scatter_size,
                                         split_contigs=args.split_contigs))
    out_files = region_lists_to_scatter_files(scattered_chunks, args.prefix)
    if args.print_paths:
        print("\n".join(out_files))
