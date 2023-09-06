import argparse

def parse_args():
	parser = argparse.ArgumentParser()
	
	parser.add_argument(
		"--genotypes", type=str, required=True, help="Path to genotypes file.")
	
	parser.add_argument(
		"--ids", type=int, nargs='+', required=True, help="Indexes of individuals that have at least one alternate allele.")

	parser.add_argument(
		"--out", type=str, required=True, help="Path to out file.")
		
	parser.add_argument(
		"--region", type=str, nargs='+', required=False, help="Region to consider formatted as chromosome, position start, and position end in 1-based fully closed coordinates.")
		
	args = parser.parse_args()
	return args

def main():
	args = parse_args()
	
	# construct match line
	match = []
	for x in range(1,131):
		if x in args.ids:
			match.append(1)
		else: 
			match.append(0)
	print(match)
	print("length of match list: %s" % len(match))

	# open out file
	out = open(args.out, 'w')
	
	# set up region variables, if included
	if args.region is not None:
		target_chr = args.region[0].strip('chr')
		target_start = int(args.region[1])
		target_end = int(args.region[2])
		
		print(target_chr, target_start, target_end)
			
	# parse genotypes file per line
	genotypes_file = open(args.genotypes, 'r')
	genotypes_lines = genotypes_file.readlines()
	for i in range(len(genotypes_lines)):
		genotype_line = genotypes_lines[i]
		genotype_line = genotype_line.strip().split('\t')
		chr = genotype_line[0]
		if chr != target_chr:
			pass
		print(i)
		
		pos = int(genotype_line[1])
		print('pos: %s' % pos)
		ref = genotype_line[2]
		alt = genotype_line[3]
		gts = list(genotype_line[4:])
		gts = [int(i) for i in gts]
		print("length of gts list: %s" % len(gts))
		# set up region variables, if included
		if args.region is not None:	
			print('region is not NONE')
			print("chr: %s, target_chr: %s" % (chr, target_chr))
			print("pos: %s, target_start: %s" % (pos, target_start))
			print(gts)
			if chr == target_chr and pos >= target_start and pos <= target_end and gts == match:
				print('chr == target_chr and pos >= target_start and pos <= target_end and gts == match')
				out.write('chr{}\t{}\t{}\t{}\t{}\n'.format(chr, pos, ref, alt, target_start))
	
		else:
			print('region is NON')
			if gts == match:
				print('gts == match')
				out.write('chr{}\t{}\t{}\t{}\n'.format(chr, pos, ref, alt))
		
	out.close()

if __name__ == '__main__':
	main()
