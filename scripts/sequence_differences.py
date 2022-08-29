# Colin M. Brand, University of California San Francisco, 06/24/2022

import argparse
import pysam

def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--intervals", type=str, required=True,
		help="Path to intervals file. Input should be two fields: 1) the contig name and 2) a comma-delimited list of start coordinates (beginning at 0)")

	parser.add_argument(
		"--sample_1", type=str, required=True, help="Path to FASTA file for first sample.")

	parser.add_argument(
		"--sample_2", type=str, required=True, help="Path to FASTA file for second sample.")

	parser.add_argument(
		"--sample_1_id", type=str, required=True, help="ID for first sample.")

	parser.add_argument(
		"--sample_2_id", type=str, required=True, help="ID for second sample.")

	args = parser.parse_args()
	return args

def main():
	args = parse_args()

	with open(f'{args.intervals}', 'r') as ints, open(f'{args.sample_1_id}_{args.sample_2_id}_sequence_differences.txt', 'w') as out:
		line = [ line.strip().split('\t') for line in ints ]
		print(line)
		for x in line:
			print(x)
			chr = x[0]
			starts = x[1].split(',')
			for x, y in zip(starts[:-1],starts[2:]):
				seq1 = pysam.FastaFile(f'{args.sample_1}').fetch(chr, int(x), int(y))
				seq2 = pysam.FastaFile(f'{args.sample_2}').fetch(chr, int(x), int(y))
				count = 0
				for a, b in zip(seq1, seq2):
					if a.upper() != b.upper():
						count = count + 1
				out.write('{}\t{}\t{}\t{}\n'.format(chr, x, y, count))

if __name__ == '__main__':
    main()
