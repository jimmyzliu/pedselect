# pedselect

# select most informative individuals to sequence from a pedigree to maximise coverage (and power) of phenotyped individuals
# single input file: "IID FATHER MOTHER SEX GENO DNA PHENO"
# defult - consider six generations

def make_dicts(ped_file):
	# make dictionary of ind:[father,mother] and ind:[child1,child2...]
	ped_dict = {i.split()[0]:[i.split()[1],i.split()[2]] for i in open(ped_file,'r')}
	parent_dict = {i:ped_dict[i] for i in ped_dict if not "0" in ped_dict[i]}
	spouse_dict = {}
	children_dict = {}
	for i in ped_dict:
		[father,mother] = ped_dict[i]
		for parent in [father,mother]:
			if parent != "0":
				if parent in children_dict:
					children_dict[parent].append(i)
				else:
					children_dict[parent] = [i]
		if father in spouse_dict:
			spouse_dict[father].append(mother)
		else:
			spouse_dict[father] = [mother]
		if mother in spouse_dict:
			spouse_dict[mother].append(father)
		else:
			spouse_dict[mother] = [father]
	for i in spouse_dict:
		spouse_dict[i] = list(set(spouse_dict[i]))
	return [parent_dict,children_dict,spouse_dict]

# get only parents and children
def get_first(ind,parent_dict,children_dict):
	rels = []
	if ind in parent_dict:
		rels = rels + parent_dict[ind]
	if ind in children_dict:
		rels = rels + children_dict[ind]
	return rels

def get_score(ind,rounds,parent_dict,children_dict,geno_inds):
	# go down s generations - children of children etc.
	geno_rels_dict = {i:[] for i in rounds}
	all_rels_dict = {i:[] for i in rounds}
	all_rels_dict[0] = [ind]
	seen_rels = []
	inds = [ind]
	for i in rounds:
		new_inds = []
		for j in inds:
			if j in geno_inds:
				continue
			if j in children_dict:
				children = children_dict[j]
				for k in children:
					new_inds.append(k)
					all_rels_dict[i].append(k)
					if k in geno_inds:
						geno_rels_dict[i].append(k)
		inds = new_inds[:]
		seen_rels = seen_rels + inds
	# go up (and down) s generations
	inds = [ind]
	if ind in parent_dict:
		parents = parent_dict[ind]
		all_rels_dict[1] = all_rels_dict[1] + parents
		if parents[0] in geno_inds:
			geno_rels_dict[1].append(parents[0])
		if parents[1] in geno_inds:
			geno_rels_dict[1].append(parents[1])
		inds = parents[:]
		for i in rounds[1:]:
			new_inds = []
			for j in inds:
				if j in geno_inds:
					continue
				if j in children_dict:
					children = children_dict[j]
					for k in children:
						if k in seen_rels:
							continue
						all_rels_dict[i].append(k)
						new_inds.append(k)
						if k in geno_inds:
							geno_rels_dict[i].append(k)
				if j in parent_dict:
					parents = parent_dict[j]
					for k in parents:
						if k in seen_rels:
							continue
						all_rels_dict[i].append(k)
						new_inds.append(k)
						if k in geno_inds:
							geno_rels_dict[i].append(k)
			inds = new_inds[:]
			seen_rels = seen_rels + inds
	for i in rounds:
		all_rels_dict[i] = list(set(all_rels_dict[i]))
		geno_rels_dict[i] = list(set(geno_rels_dict[i]))
	score = 0.0
	for i in rounds:
		c = 1.0/2**(i-1)
		score = score + c * len(geno_rels_dict[i])
	return [score,geno_rels_dict,all_rels_dict]

import sys

def main():
	# default maximum number of meiosis events to consider
	mei = 2
	# number of inds to choose for sequencing
	n_sel = 100
	# pedigree file
	ped_file = "vz.txt"
	
	args = sys.argv[1:]
	for i in range(len(args)):
		if args[i] == "-m":
			mei = int(args[i+1])
		if args[i] == "-n":
			n_sel = int(args[i+1])
		if args[i] == "-p":
			ped_file = args[i+1]
	
	rounds = range(1,mei + 1)
	
	geno_inds = []
	dna_inds = []
	pheno_inds = []
	for i in open(ped_file,'r'):
		line = i.split()
		ind = line[0]
		if line[4] == "1":
			geno_inds.append(ind)
		if line[5] == "1":
			dna_inds.append(ind)
		if line[6] == "1":
			pheno_inds.append(ind)
	
	[parent_dict,children_dict,spouse_dict] = make_dicts(ped_file)
	# do one round of trimming bottom - remove inds without children and uninform from parent_dict
	#del_inds = []
	#for i in parent_dict:
	#	if not i in children_dict:
	#		if not i in geno_inds + dna_inds + pheno_inds:
	#			del_inds.append(i)
	
	#for i in del_inds:
	#	del parent_dict[i]
	
	# first round:
	# get scores of inds with pheno and dna but no geno
	# do sequentially - assume selected inds are genotyped in subsequent rounds
	
	picked = []
	pheno_no_geno = [i for i in pheno_inds if not i in geno_inds and i in dna_inds and not i in picked]
	pheno_no_geno_no_dna = [i for i in pheno_inds if not i in geno_inds and not i in dna_inds]
	dna_no_geno_no_pheno = [i for i in dna_inds if not i in pheno_inds and not i in geno_inds]
	
	for cc in range(len(pheno_no_geno)):
		pheno_no_geno = [i for i in pheno_inds if not i in geno_inds and i in dna_inds and not i in picked]
		scores = []
		for ind in pheno_no_geno + pheno_no_geno_no_dna:
			score = get_score(ind,rounds,parent_dict,children_dict,geno_inds)[0]
			scores.append(score)		
		min_score = min(scores[:len(pheno_no_geno)])
		min_score_i = [i for i in range(len(scores)) if scores[i] == min_score and i < len(pheno_no_geno)]
		best_i = min_score_i[0]
		best_score = sum(scores) - min_score
		for i in min_score_i:
			temp_ind = pheno_no_geno[i]
			temp_scores = []
			for ind in pheno_no_geno + pheno_no_geno_no_dna:
				if ind == pheno_no_geno[i]:
					continue
				temp_scores.append(get_score(ind,rounds,parent_dict,children_dict,geno_inds + picked + [temp_ind])[0])
			if sum(temp_scores) >= best_score:
				best_score = sum(temp_scores)
				best_i = i
			# print str(sum(temp_scores)) + " " + temp_ind
		print "Adding " + pheno_no_geno[best_i] + " (" + str(best_score) + ") to list"
		geno_inds.append(pheno_no_geno[best_i])
		picked.append(pheno_no_geno[best_i])
	
	# second round
	# assume picked ids are sequenced
	# select remaining inds with genotypes but not phenotypes based on ability to impute pheno_no_geno_no_dna
	
	prev_best_score = 0.0
	for cc in range(len(pheno_no_geno_no_dna)):
		best_score = 0.0
		best_ind = "NA"
		for ind in dna_no_geno_no_pheno:
			if ind in picked:
				continue
			scores = [get_score(i,rounds,parent_dict,children_dict,geno_inds + picked + [ind])[0] for i in pheno_no_geno_no_dna]
			if sum(scores) >= best_score:
				best_score = sum(scores)
				best_ind = ind
		if best_score == prev_best_score:
			break
		picked.append(best_ind)
		print "Adding " + best_ind + " (" + str(best_score) + ") to list"
		prev_best_score = best_score

if __name__ == "__main__":
	main()
