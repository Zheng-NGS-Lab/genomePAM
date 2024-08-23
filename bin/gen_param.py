#!/usr/bin/env python
# Helper script to generate param YAML file with input flags
import argparse
import yaml

def get_args():
    '''Read arguments from terminal and return an argument object'''
    parser = argparse.ArgumentParser(description='Generate parameters.yml with command line.')
    parser.add_argument('-F', '--FQDIR', metavar='', required=True, help='Path to dir with FASTQ')
    parser.add_argument('-O', '--OUTDIR', metavar='', required=True, help='Path to output DIR')
    # parser.add_argument('-S', '--SAMPLEINFO', metavar='', required=True, help='Path to sample sheet')
    parser.add_argument('-y', '--yaml', metavar='', required=True, help='Path to YAML param file to write')
    parser.add_argument('--Read1Tail', metavar='', default='AGATCGGAAGAGCACACGTC', help='(Optional) Read1Tail')
    parser.add_argument('--Read2Tail', metavar='', default='AGATCGGAAGAGCGTCGTGT', help='(Optional) Read2Tail')
    parser.add_argument('--pos1', metavar='', default=11, help='(Optional) pos1')
    parser.add_argument('--pos2', metavar='', default=18, help='(Optional) pos2')
    parser.add_argument('--posR2', metavar='', default=8, help='(Optional) posR2')
    parser.add_argument('--xNs', metavar='', default='NNNNNNNNNN', help='(Optional) xNs')
    parser.add_argument('--FIXSEQ', metavar='', default='AGTGACAC', help='(Optional) FIXSEQ')
    parser.add_argument('--BWATHREADS', metavar='', default=4, help='(Optional) Threads for BWA')
    parser.add_argument('--GENOME', metavar='', default='/data2/database/reference_genome/human/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna', help='(Optional) Human reference genome')
    parser.add_argument('--HGVER', metavar='', default='hg38', help='(Optional) Human genome ver')
    args = parser.parse_args()
    return args


def args_to_yaml(args):
    # Convert args to dict
    yaml_dict = {
        'FQDIR': args.FQDIR, 
        'OUTDIR': args.OUTDIR,
        'Read1Tail': args.Read1Tail,
        'Read2Tail': args.Read2Tail,
        'pos1': args.pos1,
        'pos2': args.pos2,
        'posR2': args.posR2, 
        'xNs': args.xNs,
        'FIXSEQ': args.FIXSEQ,
        'BWATHREADS': args.BWATHREADS,
        'GENOME': args.GENOME,
        'HGVER': args.HGVER        
    }
    
    # Dump yaml_dict to YAML param file
    with open(args.yaml, 'w') as file:
        yaml.dump(yaml_dict, file)


if __name__ == "__main__":
    args = get_args()
    print("[LOG] Writing arguments to YAML...")
    args_to_yaml(args)
    print("[LOG] DONE!")
