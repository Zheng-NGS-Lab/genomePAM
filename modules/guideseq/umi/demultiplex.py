from __future__ import print_function
import os
import gzip
import argparse
import time
import logging

__author__ = 'Martin Aryee'

logger = logging.getLogger('root')


def fq(file):
    if file.endswith('.gz'):
        fastq = gzip.open(file, 'rt')
    else:
        fastq = open(file, 'r')
    with fastq as f:
        while True:
            l1 = f.readline()
            if not l1:
                break
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            yield [l1, l2, l3, l4]


def get_sample_id(i1, i2, sample_names):
    seq1 = i1[1]
    seq2 = i2[1]
    sample_barcode = seq1[1:8] + seq2[1:8]
    if sample_barcode in sample_names:
        return sample_names[sample_barcode]
    else:
        return sample_barcode


def demultiplex(read1, read2, index1, index2, sample_barcodes, out_dir, min_reads=10000):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if type(sample_barcodes) != dict:
        sample_names = {}
        if not sample_barcodes==None:
            with open(sample_barcodes, 'r') as barcode_file:
                for line in barcode_file:
                    fields = line.strip().split('\t')
                    if len(fields)==2:
                        sampleid, barcode = fields
                        sample_names[barcode] = sampleid
    else:
        sample_names = sample_barcodes

    outfiles_r1 = {}
    outfiles_r2 = {}
    outfiles_i1 = {}
    outfiles_i2 = {}

    total_count = 0
    count = {}
    buffer_r1 = {}
    buffer_r2 = {}
    buffer_i1 = {}
    buffer_i2 = {}

    try:
        start = time.time()
        for r1,r2,i1,i2 in zip(fq(read1), fq(read2), fq(index1), fq(index2)):
            total_count += 1
            if total_count % 1000000 == 0:
                logger.info("Processed %d reads in %.1f minutes.", total_count, (time.time()-start)/60)
            sample_id = get_sample_id(i1, i2, sample_names)

            # Increment read count and create output buffers if this is a new sample barcode
            if sample_id not in count:
                count[sample_id] = 0
                buffer_r1[sample_id] = []
                buffer_r2[sample_id] = []
                buffer_i1[sample_id] = []
                buffer_i2[sample_id] = []
            count[sample_id] += 1

            # Write the reads to file or to sample-specific output buffers if we haven't
            # reached min_reads
            # The buffers are used to avoid writing fastqs for sample barcodes with very few reads.
            if count[sample_id] < min_reads:
                buffer_r1[sample_id].append(r1)
                buffer_r2[sample_id].append(r2)
                buffer_i1[sample_id].append(i1)
                buffer_i2[sample_id].append(i2)
            elif count[sample_id] == min_reads:
                outfiles_r1[sample_id] = open(os.path.join(out_dir, '%s.r1.fastq' % sample_id), 'w')
                outfiles_r2[sample_id] = open(os.path.join(out_dir, '%s.r2.fastq' % sample_id), 'w')
                outfiles_i1[sample_id] = open(os.path.join(out_dir, '%s.i1.fastq' % sample_id), 'w')
                outfiles_i2[sample_id] = open(os.path.join(out_dir, '%s.i2.fastq' % sample_id), 'w')
                # Spill the buffers to sample-specific fastqs
                for record in buffer_r1[sample_id] + r1:
                    outfiles_r1[sample_id].write(''.join(record))
                for record in buffer_r2[sample_id] + r2:
                    outfiles_r2[sample_id].write(''.join(record))
                for record in buffer_i1[sample_id] + i1:
                    outfiles_i1[sample_id].write(''.join(record))
                for record in buffer_i2[sample_id] + i2:
                    outfiles_i2[sample_id].write(''.join(record))
                del buffer_r1[sample_id]
                del buffer_r2[sample_id]
                del buffer_i1[sample_id]
                del buffer_i2[sample_id]
            else:
                for line in r1:
                    print (line, file=outfiles_r1[sample_id], end="")
                for line in r2:
                    print (line, file=outfiles_r2[sample_id], end="")
                for line in i1:
                    print (line, file=outfiles_i1[sample_id], end="")
                for line in i2:
                    print (line, file=outfiles_i2[sample_id], end="")
    finally:
        # Close per-sample output files
        for sample_id in outfiles_r1:
            outfiles_r1[sample_id].close()
            outfiles_r2[sample_id].close()
            outfiles_i1[sample_id].close()
            outfiles_i2[sample_id].close()

    # Write remaining buffered reads to a single fastq.
    # (These reads correspond to barcodes that were seen less than min_reads times)
    with open(os.path.join(out_dir, 'undetermined.r1.fastq'), 'w') as undetermined_r1, \
         open(os.path.join(out_dir, 'undetermined.r2.fastq'), 'w') as undetermined_r2, \
         open(os.path.join(out_dir, 'undetermined.i1.fastq'), 'w') as undetermined_i1, \
         open(os.path.join(out_dir, 'undetermined.i2.fastq'), 'w') as undetermined_i2:
        for sample_id in buffer_r1.keys():
            for record in buffer_r1[sample_id]:
                undetermined_r1.write(''.join(record))
            for record in buffer_r2[sample_id]:
                undetermined_r2.write(''.join(record))
            for record in buffer_i1[sample_id]:
                undetermined_i1.write(''.join(record))
            for record in buffer_i2[sample_id]:
                undetermined_i2.write(''.join(record))

    num_fastqs = len([v for k,v in count.items() if v>=min_reads])
    logger.info('Wrote FASTQs for the %d sample barcodes out of %d with at least %d reads.', num_fastqs, len(count), min_reads)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read1', required=True)
    parser.add_argument('--read2', required=True)
    parser.add_argument('--index1', required=True)
    parser.add_argument('--index2', required=True)
    parser.add_argument('--min_reads', type=int, default=10000)
    parser.add_argument('--sample_barcodes')
    parser.add_argument('--out_dir', default='.')
    args = vars(parser.parse_args())

    demultiplex(args['read1'], args['read2'], args['index1'], args['index2'], args['sample_barcodes'], args['out_dir'], min_reads=args['min_reads'])

if __name__ == '__main__':
    main()
