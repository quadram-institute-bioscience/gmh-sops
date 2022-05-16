#!/usr/bin/env python
'''
 Script from maxibor/metaphlan-nf
 Homepage: https://github.com/maxibor/metaphlan-nf
 Author: Maxime Borry <borry@shh.mpg.de>
'''

import argparse
import pandas as pd


def _get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
        prog='metaphlan_parse',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Parsing metaphlan')
    parser.add_argument('mpReport', help="path to metaphlan output file")
    parser.add_argument(
        '-n',
        dest='nb_reads',
        default=1,
        help='Number of reads in original sample'
    )
    parser.add_argument(
        '-o',
        dest="output",
        default=None,
        help="Output file. Default = <basename>.metaphlan_parsed.csv")

    args = parser.parse_args()
    nb_reads = int(args.nb_reads)
    infile = args.mpReport
    outfile = args.output

    return(infile, nb_reads, outfile)


def _get_basename(file_name):
    if ("/") in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    return(basename)


def parse_metaphlan(infile, nb_reads, basename, outfile):
    '''
    INPUT:
        infile (str): path to metaphlan output file
        nb_reads (int): Number of reads in original sample
        basename (str): Sample basename
        outfile (str): path to output file
    OUTPUT:
        pandas Dataframe
    '''
    df = pd.read_csv(infile, sep="\t", skiprows=3)
    df.rename(columns={'#clade_name': 'Taxon',
                       'relative_abundance': basename}, inplace=True)
    df.drop(['NCBI_tax_id', 'additional_species'], axis=1, inplace=True)
    df[basename] = (df[basename]/100*nb_reads).astype(int)

    df.to_csv(outfile)


if __name__ == '__main__':
    INFILE, NB_READS, outfile = _get_args()

    basename = _get_basename(INFILE)
    if not outfile:
        outfile = basename+".metaphlan_parsed.csv"

    parse_metaphlan(
        infile=INFILE, nb_reads=NB_READS, basename=basename, outfile=outfile)