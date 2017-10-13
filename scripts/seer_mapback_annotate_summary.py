import argparse as arg

parser = arg.ArgumentParser(description='Get summary from *mapback_annotate.tsv files (by LSB)', usage = '%(prog)s [options]')
parser.add_argument("-i", "--input", help="*_mapback_annotate.tsv file", required=True)
arg = parser.parse_args()

kmerinfo = {}
geneinfo = {}
productinfo = {}
with open(arg.input, 'r') as infile:
	header = infile.readline()
	for line in infile:
		linesplit = line.rstrip().split('\t')
		kmer = linesplit[0]
		if 'NNN' not in kmer:
			maf = linesplit[1]
			wald = linesplit[2]
			score = linesplit[3]
			if len(linesplit) <= 8:
				gene = ''
				annot = ''
			else:
				gene = linesplit[11]
				annot = linesplit[15]
			l = [kmer, maf, wald, score]
			l = '\t'.join(l)
			if kmer not in kmerinfo:
				kmerinfo[kmer] = l
				geneinfo[kmer] = [gene]
				productinfo[kmer] = [annot]
			else:
				geneinfo[kmer].append(gene)
				productinfo[kmer].append(annot)

outfile = arg.input.replace('.tsv', '_summary.tsv')
with open(outfile, 'w') as out:
	header = ['kmer', 'maf', 'p_adj', 'minlog10(p_adj)', 'gene', 'annot']
	outline = '\t'.join(header)
	out.write(outline+'\n')
	for i in kmerinfo:
		info = kmerinfo[i]
		genes = set(geneinfo[i])
		annot = set(productinfo[i])
		if len(genes)>1:
			genes = 'multiple'
		else:
			genes = list(genes)[0]
		if len(annot)>1:
			annot = 'multiple'
		else:
			annot = list(annot)[0]
		outline = info+'\t'+genes+'\t'+annot+'\n'
		out.write(outline)


