==========
Changelog
==========


version 0.2.0
---------------------------
+ The ``--input`` has now been abolished in favor of a positional input
  argument.
+ The ``--prefix`` flag now has ``scatter-`` as default value.
+ Added a ``--print-paths`` or ``-P`` flag that prints the paths of the output
  files to stdout.
+ Added support for fasta index (``.fai``) files.
+ Added a ``--split-contigs`` flag to allow contigs to be split up over
  multiple files.
+ Added new ``scatter-regions`` command that is a replacement for
  `biopet-scatterregions <https://github.com/biopet/scatterregions>`_. Unlike
  biopet-scatterregions, ``scatter-regions`` works properly with scatter sizes
  larger than the Java's integer limit (2147483647). Scatter-regions invokes
  the same algorithms as ``chunked-scatter`` but with defaults that make more
  sense for a GATK scattering use case.
+ The code has been refactored for better maintainability. Python 3.5 is no
  longer supported.

version 0.1.0
---------------------------
+ First release of chunked scatter