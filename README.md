# HiTAC Assembler
HiTAC Assembler, a tool to process HiTAC-seq data

## Installing HiTAC Assembler

HiTAC Assembler is a tool for processing HiTAC-seq data in FASTQ format, which can slaso be optionall compressed by gzip.

Additional tools required for running HiTAC Assembler include:

* BWA
* SAMtools
* SPAdes

## Usage

Usage infor is as follows:

```
Usage:
HiTAC.py -t <type> -b <barcode> [options] -i <filename> -I <filename>

Input:
-t	Input data type (string). Type is one of:
	  dna - DNA fragments
	  plasmid - gene or fragment of interest in the plasmid backbone
-i <filename>	read1 input file
-I <filename>	read2 input file
-b <filename>	tn5 barcode list

Optional arguments:
-r <filename>	reference sequence filename for plasmid data
-k <int>	comma-separated list of k-mer sizes (must be odd and less than 128) [default: 91]
-T <int>	number of threads to run [default: 1]
-h	show this help message and exit
-v	show program's version number and exit
```

## HiTAC Assembler

### For plasmid data:

```bash
HiTAC.py -t plasmid -b tn5_barcode.txt -i read1.fq -I read2.fq -r reference.fa
```

### For DNA fragment data:

```bash
HiTAC.py -t dna -b tn5_barcode.txt -i read1.fq -I read2.fq
```
