import math
import argparse as arg

parser = arg.ArgumentParser(description='Convert SAM file to Manhattan plot readable by Phandango', usage = '%(prog)s [options]')
parser.add_argument("-s", "--sam", help="SAM file, resulting from mapping significant kmers to a reference", required=True)
parser.add_argument("-k", "--kmers", help="Kmers file, output from SEER", required=True)
parser.add_argument("-m", "--map", help="Mapping quality", required=False, default=0)
parser.add_argument("-o", "--outprefix", help="Output file prefix", required=True)

arg = parser.parse_args()

mappedKmers = {}
with open(arg.sam, 'r') as SAM:
	for line in SAM:
		if not line.startswith('@'):
			splitline = line.split('\t')
			pos = splitline[3]
			kmer = splitline[9]
			if splitline[2] is not '*' and int(splitline[4]) >= arg.map:
				# Get position, check if reverse complemented #
				mult = int(splitline[1]) % 16
				if mult == 0:					# if it's 0, it is multiple of 16: mapped on reverse strand
					start = int(pos)-len(kmer)+1
					end = int(pos)
				else:
					start = int(pos)
					end = int(pos)+len(kmer)-1
				# Check if multi hit, in that case the kmer may already be in the dictionary
				if kmer in mappedKmers.keys():
					mappedKmers[kmer] = mappedKmers[kmer]+"|"+str(start)+".."+str(end)
				else:
					mappedKmers[kmer] = str(start)+".."+str(end)

seerKmers = {}
with open(arg.kmers, 'r') as KMERS:
	for line in KMERS:
		kmerline = line.split('\t')
		if kmerline[0] in mappedKmers.keys():
			adj = float(kmerline[4]) # Gets LTR as adjusted p-value
			if adj==0:
				log_p = 386 # Exponent limit of a double
			elif adj>0:
				log_p = -math.log(adj)/math.log(10)
			seerKmers[kmerline[0]] = log_p

# Print output
outfile = arg.outprefix+'.phandango.plot'
with open(outfile, 'w') as out:
	for i in mappedKmers.keys():
		kmer = i
		outline = ""
		if kmer in seerKmers.keys():
			if '|' in mappedKmers[kmer]:
				curpos = mappedKmers[kmer].split('|')
				seerk = str(seerKmers[kmer])
				for cur in curpos:
					outlist = ['26', kmer, cur, seerk, '0']
					outline = '\t'.join(outlist)
					out.write(outline+'\n')
			else:
				outlist = ['26', kmer, mappedKmers[kmer], str(seerKmers[kmer]), '0']
				outline = '\t'.join(outlist)
				out.write(outline+'\n')



