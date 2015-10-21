import sys

def NX(x,lengths):
	lengths.sort(reverse=True)
	N100 = sum(lengths)
	stop = N100/ (100/float(x))
	curr_sum = 0
	for length in lengths:
		curr_sum += length
		if curr_sum >= stop:
			return length

def decide_scaffolding_procedure(Scaffolds,small_scaffolds, param):
	lengths_long = map(lambda x: Scaffolds[x].s_length, Scaffolds)
	lengths_short = map(lambda x: small_scaffolds[x].s_length, small_scaffolds)
	lengths = lengths_long + lengths_short
	N60 = NX(60,lengths)
	if N60 < param.ins_size_threshold and param.extend_paths:
		param.no_score = True
		param.contig_threshold = N60 - 1
		sys.stdout.write('Too fragmented assembly (given the insert size of the library)\
			for statistical scoring. Proceding with path finding.')
	else:
		pass
