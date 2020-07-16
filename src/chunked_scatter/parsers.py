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

import os
from typing import Generator, NamedTuple, Optional, Union

class BedRegion(NamedTuple):
    """A class that contains a region described as in the BED file format."""
    contig: str
    start: int
    end: int

    def __str__(self):
        return f"{self.contig}\t{self.start}\t{self.end}"


def dict_file_to_regions(in_file: Union[str, os.PathLike]
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


def fai_file_to_regions(in_file: Union[str, os.PathLike]
                        ) -> Generator[BedRegion, None, None]:
    # faidx format described here: https://www.htslib.org/doc/faidx.html
    with open(in_file, "rt") as in_file_h:
        for line in in_file_h:
            # faidx has name, length, offset, linebases, linewidth columns. And
            # optionally a qualoffset. By using maxsplit=2, we catch name and
            # length, while not attempting to distinguish whether we have 5 or
            # 6 columns.
            name, length, _ = line.strip().split(maxsplit=2)
            yield BedRegion(name, 0, int(length))


def file_to_regions(in_file: Union[str, os.PathLike]):
    base, extension = os.path.splitext(in_file)
    if extension == ".bed":
        return bed_file_to_regions(in_file)
    elif extension == ".dict":
        return dict_file_to_regions(in_file)
    elif extension == ".fai":
        return fai_file_to_regions(in_file)
    else:
        raise NotImplementedError(
            "Only files with .bed, .fai or .dict extensions are supported.")
