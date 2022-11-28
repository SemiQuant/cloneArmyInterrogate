# cloneArmyInterrogate

This program identifies haplotypes from Illumina paired end sequencing of a specific amplicon.

## **Usage**

python3 process_amplicons.sh "path to read directory" "path to reference bwa indexed fasta" "number of threads, def=4"

## **Required software**

-   [python3](https://www.python.org/downloads/)

    -   [pysam](https://pysam.readthedocs.io/en/latest/api.html)

    -   [pandas](https://pandas.pydata.org/getting_started.html)

    -   operator

    -   collections

-   [bwa](https://bio-bwa.sourceforge.net/)

-   [samtools](http://www.htslib.org/)
