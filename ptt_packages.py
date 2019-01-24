#!/usr/bin/env python

import os, sys
from Bio import SeqIO


def input_barcode(tn5, f, file):
    barcode = open(file, 'r')
    os.system('mkdir split_data')
    for line in barcode:
        tn5[line.split()[1]+line.split()[2]] = line.split()[0]
        f['%s%s_1' % (line.split()[1], line.split()[2])] = open('split_data/%s_1.fq' % line.split()[0], 'w')
        f['%s%s_2' % (line.split()[1], line.split()[2])] = open('split_data/%s_2.fq' % line.split()[0], 'w')
    barcode.close()


def split_data(read1, read2, tn5, f):
    print('\n\nSpliting data started.')
    import gzip
    if '.gz' in read1:
        file1 = gzip.open(read1, 'r')
    else:
        file1 = open(read1,'r')
    if '.gz' in read2:
        file2 = gzip.open(read2, 'r')
    else:
        file2 = open(read2, 'r')
    barcode_num = 0
    for r1name in file1:
        try:
            r1name = bytes.decode(r1name)
        except TypeError:
            pass
        line = next(file1)
        try:
            r1seq = bytes.decode(line)
        except TypeError:
            r1seq = line
        line = next(file1)
        try:
            r1_line3 = bytes.decode(line)
        except TypeError:
            r1_line3 = line
        line = next(file1)
        try:
            r1quality = bytes.decode(line)
        except TypeError:
            r1quality = line
        line = next(file2)
        try:
            r2name = bytes.decode(line)  # name
        except TypeError:
            r2name = line
        line = next(file2)
        try:
            r2seq = bytes.decode(line)  # seq
        except TypeError:
            r2seq = line
        line = next(file2)
        try:
            r2_line3 = bytes.decode(line)  # +
        except TypeError:
            r2_line3 = line
        line = next(file2)
        try:
            r2quality = bytes.decode(line)  # quality
        except TypeError:
            r2quality = line
        barcode_num += 1

        if tn5.__contains__(r1seq[:5] + r2seq[:5]):
            f['%s%s_1' % (r1seq[:5], r2seq[:5])].write(r1name + r1seq[24:] + '+\n' + r1quality[24:])
            f['%s%s_2' % (r1seq[:5], r2seq[:5])].write(r2name + r2seq[24:] + '+\n' + r2quality[24:])

    print(barcode_num)


def bwa(read1, read2, reference, threads):
    print('\n\nRunning BWA.')
    for i in ('.amb', '.ann', '.bwt', '.pac', '.sa'):
        if not os.path.exists(reference + i):
            os.system(f'bwa index {reference}')
            break
    os.system(f'bwa mem -t {threads} {reference} {read1} {read2} | samtools view -@ {threads} -S -f4 -f8 -b - > filter.bam')
    os.system(f'samtools fastq -1 filter.R1.fastq -2 filter.R2.fastq -n filter.bam')


def trim_data():
    print('\n\nRunning Trimmomatic.')
    g = (chr(i + 64) + str(j) for i in range(1, 9) for j in range(1, 13))
    os.system('mkdir trim_data')
    for i in g:
        if (os.path.exists(f'split_data/{i}_1.fq')) and (os.path.exists(f'split_data/{i}_2.fq')):
            os.system(f'java -jar {sys.path[0]}/Trimmomatic-0.38/trimmomatic-0.38.jar \
                      PE -phred33 split_data/{i}_1.fq split_data/{i}_2.fq \
                      {i}_1.paired.fq {i}_1.unpaired.fq {i}_2.paired.fq {i}_2.unpaired.fq \
                      ILLUMINACLIP:/NAS7/home/gaoxiang/softwares/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:30:10 \
                      SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:112')
            os.system('mv *paired.fq trim_data')


def assemble(types, kmer):
    print('\n\nRunning SPAdes.')
    g = (chr(i + 64) + str(j) for i in range(1, 9) for j in range(1, 13))
    if types == 'dna':
        for i in g:
            if (os.path.exists(f'trim_data/{i}_1.paired.fq')) and (os.path.exists(f'trim_data/{i}_2.paired.fq')):
                os.system(f'spades.py -1 trim_data/{i}_1.paired.fq -2 trim_data/{i}_2.paired.fq --careful -k {kmer} -o spades/{i}')
    if types == 'plasmid':
        for i in g:
            if (os.path.exists(f'split_data/{i}_1.fq')) and (os.path.exists(f'split_data/{i}_2.fq')):
                os.system(f'spades.py -1 split_data/{i}_1.fq -2 split_data/{i}_2.fq --careful -k {kmer} -o spades/{i}')


def select_contig():
    import re
    g = (chr(i + 64) + str(j) for i in range(1, 9) for j in range(1, 13))
    os.system('mkdir ptt_output')
    output_file = open('ptt_output/ptt_output.fa', 'w')
    for i in g:
        if os.path.exists(f'spades/{i}/scaffolds.fasta'):
            os.system(f'cp spades/{i}/scaffolds.fasta ptt_output/{i}.fa')
            max_score = 0
            max_name = ''
            max_seq = ''
            for seq in SeqIO.parse(f'ptt_output/{i}.fa', 'fasta'):
                name = seq.id
                seq = str(seq.seq)
                num = re.findall(r'\d+\.?\d*', name)
                if int(num[1]) * float(num[2]) > max_score:
                    max_score = int(num[1]) * float(num[2])
                    max_name = name
                    max_seq = seq
            if max_score != 0:
                output_file.write(f'>{i} {max_name}\n{max_seq}\n')
    output_file.close()
    result_path = os.popen('pwd').readline()[:-1]
    print(f'\n\nPTT-seq result are at: {result_path}/ptt_output/')
