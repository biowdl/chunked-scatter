# chunked_scatter

This tool takes a bed file or sequence dictionary as input and divides the
contigs/chromosomes into overlapping chunks of a given size. These chunks will
then be placed in new bed files, one chromosomes per file. Small chromosomes
will be put together to avoid the creation of thousands of files.

## Installation
Install from github:
- Clone the repository: `git clone https://github.com/biowdl/chunked-scatter.git`
- Enter the repository: `cd chunked-scatter`
- Install using pip: `pip install .`

## Usage
```
chunked-scatter -p output_prefix -i input.bed
```
The input is expected to end in `.bed` or `.dict`!

| option | arguments | definition |
|-|-|-|
| -c | a number | The size of the chunks. |
| -o | a number | The size of the overlap. |
| -m | a number | The minimum number of bases to be put in a single output file, before a new scatter will be made. |

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
chunked-scatter -p /data/scatter_ -i /data/regions.bed -m 1000 -c 5000
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
chunked-scatter -p /data/scatter_ -i /data/regions.bed
```

Will produce the following output file at `/data/scatter_0.bed`:
```
chr1	0	1000000
chr1	999850	2000000
chr1	1999850	3000000
chr2	0	500000
```
