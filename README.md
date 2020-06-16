# chunked-scatter and scatter-regions

The `chunked-scatter` tool takes a bed file or sequence dictionary as input 
and divides the
contigs/chromosomes into overlapping chunks of a given size. These chunks will
then be placed in new bed files, one chromosomes per file. Small chromosomes
will be put together to avoid the creation of thousands of files.

The `scatter-regions` tool works in a similar way but with defaults and flags
tuned towards creating genome scatters for GATK tools.

## Installation
- Install using pip: `pip install chunked-scatter`
- Install using conda: `conda install chunked-scatter`
    - This requires [conda with a bioconda channel](
    http://bioconda.github.io/user/install.html#).

## Usage

### chunked-scatter
```
usage: chunked-scatter [-h] [-p PREFIX] [-S] [-P] [-c SIZE]
                       [-m MINIMUM_BP_PER_FILE] [-o OVERLAP]
                       INPUT

Given a sequence dict, fasta index or a bed file, scatter over the defined
contigs/regions. Each contig/region will be split into multiple overlapping
regions, which will be written to a new bed file. Each contig will be placed
in a new file, unless the length of the contigs/regions doesn't exceed a given
number.

positional arguments:
  INPUT                 The input file, either a bed file or a sequence dict.
                        Which format is used is detected by the extension:
                        '.bed', '.fai' or '.dict'. This option is mandatory.

optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        The prefix of the ouput files. Output will be named
                        like: <PREFIX><N>.bed, in which N is an incrementing
                        number. Default 'scatter-'.
  -S, --split-contigs   If set, contigs are allowed to be split up over
                        multiple files.
  -P, --print-paths     If set prints paths of the output files to STDOUT.
                        This makes the program usable in scripts and
                        worfklows.
  -c SIZE, --chunk-size SIZE
                        The size of the chunks. The first chunk in a region or
                        contig will be exactly length SIZE, subsequent chunks
                        will SIZE + OVERLAP and the final chunk may be
                        anywhere from 0.5 to 1.5 times SIZE plus overlap. If a
                        region (or contig) is smaller than SIZE the original
                        regions will be returned. Defaults to 1e6
  -m MINIMUM_BP_PER_FILE, --minimum-bp-per-file MINIMUM_BP_PER_FILE
                        The minimum number of bases represented within a
                        single output bed file. If an input contig or region
                        is smaller than this MINIMUM_BP_PER_FILE, then the
                        next contigs/regions will be placed in the same file
                        untill this minimum is met. Defaults to 45e6.
  -o OVERLAP, --overlap OVERLAP
                        The number of bases which each chunk should overlap
                        with the preceding one. Defaults to 150.
```

### scatter-regions
```
usage: scatter-regions [-h] [-p PREFIX] [-S] [-P] [-s SCATTER_SIZE] INPUT

Given a sequence dict, fasta index or a bed file, scatter over the defined
contigs/regions. Creates a bed file where the contigs add up approximately to
the given scatter size.

positional arguments:
  INPUT                 The input file, either a bed file or a sequence dict.
                        Which format is used is detected by the extension:
                        '.bed', '.fai' or '.dict'. This option is mandatory.

optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        The prefix of the ouput files. Output will be named
                        like: <PREFIX><N>.bed, in which N is an incrementing
                        number. Default 'scatter-'.
  -S, --split-contigs   If set, contigs are allowed to be split up over
                        multiple files.
  -P, --print-paths     If set prints paths of the output files to STDOUT.
                        This makes the program usable in scripts and
                        worfklows.
  -s SCATTER_SIZE, --scatter-size SCATTER_SIZE
                        How large the regions over which to scatter should be.
                        Default: 1000000000.

```

## Examples
### bed file
Given a bed file located at `/data/regions.bed`:
```
chr1	100	1000
chr1	2000	16000
chr2	5000	10000
```

The command:
```
chunked-scatter -p /data/scatter_ -m 1000 -c 5000 /data/regions.bed
```

Will produce the following two output files:
- `/data/scatter_0.bed`:
  ```
  chr1	100	1000
  chr1	2000	7000
  chr1	6850	12000
  chr1	11850	16000
  ```
- `/data/scatter_1.bed`:
  ```
  chr2	5000	10000
  ```

### dict file
Given a dict file located at `/data/ref.dict`:
```
@SQ	SN:chr1	LN:3000000
@SQ SN:chr2 LN:500000
```

The command:
```
chunked-scatter -p /data/scatter_ /data/regions.bed
```

Will produce the following output file at `/data/scatter_0.bed`:
```
chr1	0	1000000
chr1	999850	2000000
chr1	1999850	3000000
chr2	0	500000
```
