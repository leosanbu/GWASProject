import time
import math
import random
import argparse as arg
from Bio import SeqIO

parser = arg.ArgumentParser(description='Locate and annotate kmers from SEER (by LSB)', usage = '%(prog)s [options]')
parser.add_argument("-k", "--kmers", help="Kmers file, filtered output from SEER (required)", required=True)
parser.add_argument("-s", "--assembly", help="File containing two columns: strain and location of assembly (required)", required=True)
parser.add_argument("-a", "--annot", help="File containing two columns: strain and location of GFF files (optional)", required=False)
parser.add_argument("-o", "--outprefix", help="Output file prefix (required)", required=True)
parser.add_argument("-v", "--pvalue", help="Use Wald (continuous data, 0) or LTR (binary data, 1) p-values", required=False, default=0)
parser.add_argument("-f", "--filter_pval", help="Filter out kmers with an adjusted p-value (from the one selected with -v) above this threshold (default = 1, do not filter)", required=False, default=1)
parser.add_argument("-m", "--max", help="Maximum number of gff files to look at for annotating a kmer (optional, default = 5)", required=False, default=5)
parser.add_argument("-l", "--list", help="File containing a list of strain names with good annotations to always keep (must be same name as the first column in -s and -a) (optional)", required=False)
parser.add_argument("-c", "--column", help="Column number where the list of strains containing a kmer starts (optional, default = 8)", required=False, default=8)
parser.add_argument("-z", "--summary", help="Create summary of final output (optional) (default=True)", required=False, default=True)
parser.add_argument("-p", "--phandango", help="Create a Manhattan plot input file for Phandango from this strain", required=False)
parser.add_argument("-t", "--time", help="Print timing on screen every this number of files (optional, default = 100)", required=False, default=100)
arg = parser.parse_args()

starttime = time.time()
prev = 0

## FUNCTIONS ##

def readLocations(filename):
	dicLoc = {}
	with open(filename, 'r') as LOC:
		for line in LOC:
			linesplit = line.rstrip().split('\t')
			dicLoc[linesplit[0]] = linesplit[1]
	return dicLoc

def parseGFF(gff):
	gffDic = {}
	gffCoords = {}
	with open(gff, 'rU') as gffile:
		for entry in gffile:
			if not entry.startswith('##') and '\t' in entry:
				entrysplit = entry.rstrip().split('\t')
				contig = entrysplit[0]
				ini = int(entrysplit[3])
				end = int(entrysplit[4])
				strand = entrysplit[6]
				annot = strand+';'+entrysplit[8]
				if contig in gffCoords:
					gffCoords[contig][ini] = end
					gffDic[contig][ini] = annot
				else:
					gffCoords[contig] = {}
					gffCoords[contig][ini] = end
					gffDic[contig] = {}
					gffDic[contig][ini] = annot
	return [gffDic, gffCoords]

def checkGene(contig, gene, gffDic, gffCoords, closest_match):
	positions = list(gffDic[contig].keys())
	gene_annot = gffDic[contig][positions[gene]] 
	gene_strand = gene_annot[0]
	gene_start = positions[gene]
	gene_end = gffCoords[contig][positions[gene]]
	# upstream or inside the closest gene?
	where = ""
	if gene_strand == '+' and closest_match<0:
		where = "upstream"
	elif gene_strand == '+' and closest_match>=0:
		suma = gene_start+closest_match
		if suma > gene_end:
			where = "downstream"
		else:
			where = "coding"
	elif gene_strand == '-' and closest_match<0:
		suma = gene_end+closest_match
		if suma < gene_start:
			where = "downstream"
		else:
			where = "coding"
	elif gene_strand == '-' and closest_match>=0:
		where = "upstream"
	return where

def annotateKmer(start, contig, gffDic, gffCoords):
	positive = []
	negative = []
	allpos = []
	for positionStart in gffDic[contig]:
		positionEnd = gffCoords[contig][positionStart]
		positionStrand = gffDic[contig][positionStart][0]
		originalStart = positionStart
		if positionStrand is '-':
			positionStart = positionEnd
			positionEnd = originalStart
		calc = start-positionStart
		allpos.append(calc)
		if calc<0:
			negative.append(calc)
		elif calc>=0:
			positive.append(calc)
	gene = 0
	minpositive = 0
	maxnegative = 0
	which = 0
	closest_match = 0
	where = ''
	if len(positive)==0:
		maxnegative = max(negative)
		gene = allpos.index(maxnegative)	# all values are negative
		which = -1
		closest_match = maxnegative
		where = checkGene(contig, gene, gffDic, gffCoords, maxnegative)
	elif len(negative)==0:
		minpositive = min(positive)
		gene = allpos.index(minpositive)	 # all values are positive
		which = 1
		closest_match = minpositive
		where = checkGene(contig, gene, gffDic, gffCoords, minpositive)
	else:
		minpositive = min(positive) # how far is from start position 5'->3'
		maxnegative = max(negative) # how far is from start position 5'<-3'
		poscheck = checkGene(contig, allpos.index(minpositive), gffDic, gffCoords, minpositive)
		negcheck = checkGene(contig, allpos.index(maxnegative), gffDic, gffCoords, maxnegative)
		# which gene is closest
		if minpositive < abs(maxnegative):	# minpositive wins UNLESS maxnegative is in a gene
			if poscheck is not 'coding' and negcheck is 'coding':
				gene = allpos.index(maxnegative)
				which = -1
				where = negcheck
				closest_match = maxnegative
			else:
				gene = allpos.index(minpositive)
				which = 1
				where = poscheck
				closest_match = minpositive
		elif abs(maxnegative) < minpositive: # maxnegative wins UNLESS minpositive is in a gene
			if negcheck is not 'coding' and poscheck is 'coding':
				gene = allpos.index(minpositive)
				which = 1
				where = poscheck
				closest_match = minpositive
			else:
				gene = allpos.index(maxnegative)
				which = -1
				where = negcheck
				closest_match = maxnegative
	# get gene info
	positions = list(gffDic[contig].keys())
	gene_annot = gffDic[contig][positions[gene]]
	gene_strand = gene_annot[0]
	gene_start = positions[gene]
	gene_end = gffCoords[contig][positions[gene]]
	# split annot
	genesplit = gene_annot.split(';')
	gene_locustag = ""
	gene_gene = ""
	gene_inference = ""
	gene_product = ""
	for i in genesplit:
		if i.startswith('gene='):
			gene_gene = i.replace('gene=', '').replace("\"", "")
		elif i.startswith('locus_tag='):
			gene_locustag = i.replace('locus_tag=', '').replace("\"", "")
		elif i.startswith('product='):
			gene_product = i.replace('product=', '').replace("\"", "")
		elif i.startswith('inference'):
			gene_inference = i.replace('inference=', '').replace("\"", "")
	# prepare final output
	outannot = [str(closest_match+1), where, gene_locustag, gene_gene, str(gene_start), str(gene_end), gene_strand, gene_product, gene_inference]
	return outannot

def revComp(kmer):
	basesDic = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
	kmerlist = list(kmer)
	revkmer = list(reversed(kmerlist))
	comp = [basesDic[x] for x in revkmer]
	return ''.join(comp)

def getKmerInfo(kmer, seq, gffDic, gffCoords, forward):
	seqsplit = seq.split(kmer)
	kmerlen = len(kmer)
	first = len(seqsplit[0])
	lensplit = [first]
	for x in range(1,len(seqsplit)):
		lensplit.append(len(seqsplit[x])+first+kmerlen)
		first += len(seqsplit[x])+kmerlen
	lensplit = lensplit[:-1] # ignore last one
	pos = 0
	outputlist = []
	for chunk in lensplit:
		start = chunk
		end = chunk+kmerlen
		pos = end
		pr = [name, str(start+1), str(end), str(forward)]
		if do_annot is True:
			# Annotate #
			if name in gffCoords:
				out = annotateKmer(start, name, gffDic, gffCoords)
				if forward == 1:
					output = kmer+'\t'+kmerSign[kmer]+'\t'+'\t'.join(pr)+'\t'+'\t'.join(out)+'\n'
				else:
					output = revComp(kmer)+'\t'+kmerSign[revComp(kmer)]+'\t'+'\t'.join(pr)+'\t'+'\t'.join(out)+'\n'
			else:
				# tweak contig name 
				contig_number = name.split('.')[-1]
				contig_number = '00'+contig_number
				real_contig = ''
				contigfound = 0
				for cont in gffCoords:
					if cont.endswith(contig_number):
						real_contig = cont
						contigfound = 1
				if real_contig in gffCoords:
					out = annotateKmer(start, real_contig, gffDic, gffCoords)
					if forward == 1:
						output = kmer+'\t'+kmerSign[kmer]+'\t'+'\t'.join(pr)+'\t'+'\t'.join(out)+'\n'
					else:
						output = revComp(kmer)+'\t'+kmerSign[revComp(kmer)]+'\t'+'\t'.join(pr)+'\t'+'\t'.join(out)+'\n'
				if contigfound == 0:
					if forward == 1:
						output = kmer+'\t'+kmerSign[kmer]+'\t'+'\t'.join(pr)+'\n'
					else:
						output = revComp(kmer)+'\t'+kmerSign[revComp(kmer)]+'\t'+'\t'.join(pr)+'\n'
			outputlist.append(output)
	return outputlist

# Check is annotation required #
if arg.annot is None:
	do_annot = False
else:
	do_annot = True

# Read file locations #
# Make sure all assemblies have its corresponding annotation #
assemblyLoc = readLocations(arg.assembly)
if do_annot is True:
	annotLoc = readLocations(arg.annot)
	notFound = []
	for i in assemblyLoc:
		if i not in annotLoc:
			notFound.append(i)
	if len(notFound)>0:
		with open('annotation_not_found.txt', 'w') as NOT:
			for i in notFound:
					NOT.write(i+"\n")
		raise Exception ('Annotation not found for', notFound)

# Check if there is a list of trusted annotations #
if arg.list is None:
	trusted = False
else:
	trusted = True
	trustedList = []
	with open(arg.list, 'r') as LIST:
		for line in LIST:
			trustedList.append(line.rstrip())

max_samp = int(arg.max)

comments_field = ['bfgs-fail','nr-fail','firth-fail','large-se','bad-chisq','inv-fail','zero-ll']

# Save all kmers to search in each assembly #
allsamples = {}
kmerSign = {}
byAssembly = {}
columnIndex = arg.column-1
with open(arg.kmers, 'r') as KMERS:	
	for line in KMERS:
		linesplit = line.rstrip().split('\t')
		kmer = linesplit[0]
		maf = linesplit[1]
		p_adj = float(linesplit[3])
		if int(arg.pvalue) == 1:
			p_adj = float(linesplit[4])
		beta = linesplit[5]
		ranStrains = []
		if p_adj <= float(arg.filter_pval):
			if p_adj == 0:
				score = 386
			elif p_adj > 0:
				score = -math.log(p_adj)/math.log(10)
			kmerSign[kmer] = str(maf)+'\t'+str(p_adj)+'\t'+str(float(round(score,4)))+'\t'+str(beta)
			samples = linesplit[columnIndex:]
			cleansamples = []
			for x in samples:
				if x != 'NA':
					curx = x.split(',')
					for y in curx:
						if y not in comments_field:
							cleansamples.append(x)
			allsamples[kmer] = cleansamples
			if trusted is True:
				excSamples = []
				incSamples = []
				for s in cleansamples:
					if s not in trustedList:
						excSamples.append(s)
					else:
						incSamples.append(s)
				if len(incSamples)==max_samp:
					ranStrains = incSamples
				elif len(incSamples)>max_samp:
					ranStrains = random.sample(incSamples, max_samp)
				elif len(incSamples)<max_samp:
					left = max_samp-len(incSamples)
					leftStrains = random.sample(excSamples, left)
					ranStrains = incSamples
					for x in leftStrains:
						ranStrains.append(x)
				# Make sure the key strain for phandango (if requested) is selected here #
				if arg.phandango is not None:
					if arg.phandango in cleansamples and arg.phandango not in ranStrains:
						ranStrains[0] = arg.phandango
			else:
				if len(cleansamples)<max_samp:
					max_samp = len(cleansamples)
				ranStrains = random.sample(cleansamples, max_samp)
			# Dictionary: kmers by assembly
			for i in ranStrains:
				if i in byAssembly:
					byAssembly[i].append(kmer)
				else:
					byAssembly[i] = []
					byAssembly[i].append(kmer)

count = 0
outputDic = {}
for i in byAssembly:
	strain = i
	# add to counter and print if multiple of 10 #
	count = count+1
	checkcount = count % arg.time == 0
	if checkcount is True:
		curtime = (time.time()-starttime)-prev
		prev += curtime
		print(str(count)+' in '+str(round(curtime,2))+' seconds')
	# process
	path = assemblyLoc[strain]
	gff = annotLoc[strain]
	if strain in byAssembly:
		if do_annot is True:
			# Read gff file - just once #
			getgff = parseGFF(gff)
			gffDic = getgff[0]
			gffCoords = getgff[1]
			# Read fasta file #
		with open(path, 'rU') as fasta:
			for record in SeqIO.parse(fasta, "fasta"):
				name = record.id
				seq = str(record.seq)
				for kmer in byAssembly[strain]:
					kmerlen = len(kmer)
					kmer_rc = revComp(kmer)
					if kmer in seq:
						forward = 1
						output = getKmerInfo(kmer, seq, gffDic, gffCoords, forward)
						for o in output:
							if kmer in outputDic:
								outputDic[kmer].append(o)
							else:
								outputDic[kmer] = []
								outputDic[kmer].append(o)
					# Check if there are matches in reverse comlementary
					if kmer_rc in seq:
						forward = -1
						output = getKmerInfo(kmer_rc, seq, gffDic, gffCoords, forward)
						for o in output:
							if kmer in outputDic:
								outputDic[kmer].append(o)
							else:
								outputDic[kmer] = []
								outputDic[kmer].append(o)

# Generate raw output file #
outfile = arg.outprefix+'_mapback_annotate.tsv'
with open(outfile, 'w') as OUT:
	if do_annot is True:
		header = ['kmer', 'maf', 'p_adj_wald', 'neglog10(p_adj_wald)', 'beta', 'contig', 'kmer_start', 'kmer_end', 'contig_strand', 'distance_from_start', 'where_matches', 'locus_tag', 'gene', 'gene_start', 'gene_end', 'gene_strand', 'product', 'inference']
	else:
		header = ['kmer', 'maf', 'p_adj_wald', 'neglog10(p_adj_wald)', 'beta', 'contig', 'kmer_start', 'kmer_end', 'contig_strand']
	OUT.write('\t'.join(header)+'\n')
	for kmer in outputDic:
		for result in outputDic[kmer]:
			OUT.write(result)

# Create summary file if requested #
if arg.summary is True:		#arg.summary
	kmerinfo = {}
	for entry in outputDic:
		if 'NNN' not in entry:
			geneinfo = []
			productinfo = []
			kmerList = outputDic[entry]
			for k in kmerList:
				linesplit = k.rstrip().split('\t')
				kmer = linesplit[0]
				maf = linesplit[1]
				wald = linesplit[2]
				score = linesplit[3]
				beta = linesplit[4]
				if len(linesplit) <= 9:
					gene = ''
					annot = ''
				else:
					gene = linesplit[12]
					annot = linesplit[16]
				l = [kmer, maf, wald, score, beta]
				l = '\t'.join(l)
				geneinfo.append(gene)
				productinfo.append(annot)
			# Summarize gene and product info for this kmer #
			finalgene = ''
			finalannot = ''
			# Gene #
			geneclean = [x.split('_')[0] for x in geneinfo]
			geneclean = [x for x in geneclean if x is not '']
			if len(geneclean)>0:
				uniq = list(set(geneclean))
				if len(uniq)==1:
					finalgene = uniq[0]
				else:
					# count how many of each #
					howmany = {}
					for x in uniq:
						hmany = geneclean.count(x)
						howmany[x] = hmany
					maxvalue = sorted(howmany.values())
					#if len(maxvalue)-1 < len(maxvalue) or len(maxvalue)-1 > len(maxvalue):
					#	raise Exception('Stopped in '+entry)
					maxvalue = maxvalue[len(maxvalue)-1]
					selgenes = [x for x in howmany if howmany[x]==maxvalue]
					if len(selgenes)==1:
						# if there is an option with a maximum number #
						finalgene = selgenes[0]
					else:
						finalgene = ','.join(selgenes)
			# Annotation #
			uniq = list(set(productinfo))
			if len(uniq)==1:
				finalannot = uniq[0]
			else:
				# count how many of each #
				howmany = {}
				for x in uniq:
					hmany = productinfo.count(x)
					howmany[x] = hmany
				maxvalue = sorted(howmany.values())
				#if len(maxvalue)-1 < len(maxvalue) or len(maxvalue)-1 > len(maxvalue):
				#		raise Exception('Stopped in '+entry)
				maxvalue = maxvalue[len(maxvalue)-1]
				selannot = [x for x in howmany if howmany[x]==maxvalue]
				if len(selannot)==1:
					# if there is an option with a maximum number #
					finalannot = selannot[0]
				else:
					finalannot = ','.join(selannot)
			l = l+'\t'+finalgene+'\t'+finalannot
			kmerinfo[kmer] = l

	# Write output file #
	outfile = arg.outprefix+'_mapback_summary.tsv'
	with open(outfile, 'w') as SUMM:
		header = ['kmer', 'maf', 'p_adj_wald', 'neglog10(p_adj_wald)', 'beta', 'gene', 'annotation']
		SUMM.write('\t'.join(header)+'\n')
		for i in kmerinfo:
			SUMM.write(kmerinfo[i]+'\n')

# Create Manhattan plot file for Phandango if specified
if arg.phandango is not None:
	strain = arg.phandango
	kmers_in_strain = byAssembly[strain]
	outfile = arg.outprefix+'_'+strain+'_phandango.plot'
	with open(outfile, 'w') as out:
		for k in kmers_in_strain:
			res = outputDic[k]
			ressplit = res[0].split('\t')
			score = ressplit[3]
			k_start = ressplit[6]
			k_end = ressplit[7]
			phanline = ['26', k, k_start+'..'+k_end, score, '0']
			outline = '\t'.join(phanline)
			out.write(outline+'\n')
			#26	TTTTTTTTTTTAAC	2112792..2112805	8.940058111938043	0

endtime = time.time()-starttime
print('Script finished in', round(endtime/60,2), 'minutes')

