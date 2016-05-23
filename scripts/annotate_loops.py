# annotate_loops.py
# Usage: python annotate_loops.py <loop annotations (BED4)> <loop file>
# TODO:
# 	not getting relative position of loop around DHS site yet, just getting counts

import csv
import sys
import math
import collections
import cPickle as pickle
csv.register_dialect("textdialect", delimiter='\t')

if len(sys.argv) != 4:
	print "Usage: python annotate_loops.py <loop annotations (BED4 or pkl)> <loop file> <pickle output>"
	exit()

def lx():
	return {}
	
if '.pkl' in sys.argv[1]:
	pos_to_dhs = pickle.load(open(sys.argv[1], 'r'))
else:
	ifile = open(sys.argv[1], 'r')  # "interm/nb2/v65_cohesin_merge_allValidPairs_gt20kb.sorted_dhs.bed"
	reader = csv.reader(ifile, 'textdialect')
		
	pos_to_dhs = collections.defaultdict(lx)  # chr --> {pos (rounded to 100) --> dhs}
	chr = "abc"
	ctr = 0
	for row in reader:
		ctr += 1
		if ctr % 1000000 == 0: print "Position line {}".format(ctr)
		if row[0] != chr:
			chr = row[0]
			print chr
		pos_to_dhs[row[0]][int(row[1])/100] = '__'.join(row[3].split('__')[4:6])
		
	ifile.close()
	pickle.dump(pos_to_dhs, open(sys.argv[1][:-4] + '.pkl', 'w'))
	
ifile = open(sys.argv[2], 'r')  # "data/loops/v65_cohesin_merge_allValidPairs"
reader = csv.reader(ifile, 'textdialect')

def l():
	return 0
	
def ll():
	return collections.defaultdict(l)

dhs_to_dhs = collections.defaultdict(ll)  # dhs --> {dhs --> ct}
ctr = 0
twentyctr = 0
total_good = 0
for row in reader:
	ctr += 1
	if ctr % 1000000 == 0: print "Pair line {}".format(ctr)
	if math.fabs(int(row[5]) - int(row[2])) < 20000: continue
	twentyctr += 1
	if twentyctr % 1000000 == 0: print "Pair line (for >20kb pairs) {}".format(twentyctr)
	try: 
		dhs1 = pos_to_dhs[row[1]][int(row[2])/100]
		dhs2 = pos_to_dhs[row[4]][int(row[5])/100]
	except: 
		continue
	if dhs1 < dhs2: dhs_to_dhs[dhs1][dhs2] += 1
	else: dhs_to_dhs[dhs2][dhs1] += 1
	total_good += 1
	if total_good % 10000 == 0: print "Total good pairs: {}".format(total_good)

	
ifile.close()
pickle.dump((dict(dhs_to_dhs), total_good, twentyctr, ctr), open(sys.argv[3], 'w'))
	

	



