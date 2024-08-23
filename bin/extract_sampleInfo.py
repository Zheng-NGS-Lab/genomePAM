#!/usr/bin/env python
import argparse
import pandas as pd
import sys

def get_args():
    '''Read arguments from terminal and return an argument object'''
    parser = argparse.ArgumentParser(description='Extract sample information from CSV.')
    parser.add_argument('-s', '--sample_info', metavar='', required=True, help='Sample Info CSV used in demultiplexing')
    parser.add_argument('-l', '--lib', metavar='', required=True, help='Library/Sample ID')
    parser.add_argument('-c', '--column', metavar='', required=True, help='Column name in sample info CSV to extract')

    args = parser.parse_args()
    return args


def read_sample_info(sample_info_csv:str)->pd.DataFrame:
    sample_info_df = pd.read_csv(sample_info_csv, skiprows=1, index_col=[0])
    return sample_info_df


def main():
    # Read arguments from terminal
    args = get_args()

    # Read sample info CSV as a pd.DatFrame
    sample_info_df = read_sample_info(args.sample_info)

    # Query required info from sample_info_df
    try:
        if args.lib.startswith('Lib'):
            # Sample for iseq/nextseq
            lib_name = args.lib.split('_')[0]
            required_info = sample_info_df.loc[lib_name, args.column]
        else:
            # Sample from CRO
            lib_name = args.lib
            for Sample_ID, row in sample_info_df.iterrows():
                if lib_name.startswith(Sample_ID):
                    required_info = row[args.column]
        sys.stdout.write(required_info)
    except Exception as e:
        print(e)
        sys.stderr.write(f'Error occured when fetching:\n{e}\n')


if __name__ == "__main__":
    main()
