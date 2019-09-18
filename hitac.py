#!/usr/bin/env python3

import getopt, sys, os


def usage():
    print('Usage:\n'
          'ptt.py -t <type> -b <barcode> [options] -i <filename> -I <filename>\n\n'
          'Input:\n'
          '-t\tInput data type (string). Type is one of:\n'
          '\t  dna - DNA fragments\n'
          '\t  plasmid - gene or fragment of interest in the plasmid backbone\n'
          '-i <filename>\tread1 input file\n'
          '-I <filename>\tread2 input file\n'
          '-b <filename>\ttn5 barcode list\n\n'
          'Optional arguments:\n'
          '-r <filename>\treference sequence filename for plasmid data\n'
          '-k <int>\tcomma-separated list of k-mer sizes (must be odd and less than 128) [default: 91]\n'
          '-T <int>\tnumber of threads to run [default: 1]\n'
          '-h\tshow this help message and exit\n'
          '-v\tshow program\'s version number and exit\n')


def main():

    types, read1, read2 = '', '', ''
    threads, kmer, reference, barcode = 1, 91, '', ''

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hi:I:t:k:vr:b:T:')
        # argv = 'test.py -t plasmid -i 9R-1_R1.fq.gz -I 9R-1_R2.fq.gz -k 91 -b tn5_barcode.txt -r tn5_barcode.txt'.split()
        # opts, args = getopt.getopt(argv[1:], 'hi:I:t:k:vr:b:T:')
    except getopt.GetoptError as err:
        # print help information and exit:
        print(f'Error:{err}\n')  # will print something like 'option -a not recognized'
        usage()
        sys.exit(2)

    if len(opts) == 0:
        print('Error: Option not recognized\n')
        usage()
        sys.exit(2)

    for o, a in opts:
        if o == '-v':
            print('PTT-seq v1.0.0')
            sys.exit(0)
        elif o == '-h':
            usage()
            sys.exit()
        elif o == '-t':
            types = a
        elif o == '-i':
            read1 = a
        elif o == '-I':
            read2 = a
        elif o == '-k':
            try:
                kmer = int(a)
            except ValueError:
                print(f'KmerError: \'{a}\' is not INT')
        elif o == '-r':
            reference = a
        elif o == '-b':
            barcode = a
        elif o == '-T':
            try:
                threads = int(a)
            except ValueError:
                print(f'ThreadsError: \'{a}\' is not INT')
        else:
            assert False, 'unhandled option'

    if not (types in ('dna', 'plasmid')):
        print(f'Error: Type {types} not recognized\n')
        usage()
        sys.exit(2)
    if not os.path.exists(read1):
        print('Error: Fail to locate the read1 file\n')
        sys.exit(2)
    if not os.path.exists(read2):
        print('Error: Fail to locate the read2 file\n')
        sys.exit(2)
    if os.path.exists(barcode):
        from ptt_packages import input_barcode
        tn5, output_file = {}, {}
        input_barcode(tn5, output_file, barcode)
    else:
        print('Error: Fail to locate the barcode file\n')
        sys.exit(2)

    if types == 'dna':
        from ptt_packages import split_data, trim_data, assemble
        split_data(read1, read2, tn5, output_file)
        trim_data()
        assemble(types, kmer)
    if types == 'plasmid':
        if not os.path.exists(reference):
            print('Error: Fail to locate the reference file\n')
            sys.exit(2)
        from hitac_packages import bwa, split_data, assemble
        bwa(read1, read2, reference, threads)
        split_data('filter.R1.fastq', 'filter.R2.fastq', tn5, output_file)
        assemble(types, kmer)

    from ptt_packages import select_contig
    select_contig()

    print('\n\nHiTAC assembler finished!')


if __name__ == '__main__':
    main()

