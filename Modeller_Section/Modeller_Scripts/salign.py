# Illustrates the SALIGN multiple structure/sequence alignment

from modeller import *

import argparse

parser = argparse.ArgumentParser(description="Multiple structure/sequence alignment")
parser.add_argument('-l', '--templateList', metavar='', nargs='+', required=True, help="List of templates passed in as a bash array")
args = parser.parse_args()


def cleaning_parser_list(aString):
	aString = aString.replace('[','')
	aString = aString.replace(']','')
	aString = aString.replace(',','')
	aString = aString.replace(' ','')
	aString = aString.replace('.pdb','_A')
	return aString
	
	
def main(args):	
	log.verbose()
	env = Environ()
	env.io.atom_files_directory = ['.', '../atom_files/']
	
	working_pdb_list = []
	for pdb in args.templateList:
		working_pdb_list.append(f'{cleaning_parser_list(pdb)}')
	
	
	
	############--modified saling--############
	code_and_chain = []
	for template in working_pdb_list:
		tmp = template.split('_')
		print(tmp)
		try:
			code_and_chain.append((tmp[0],tmp[1]))
		except:
			code_and_chain.append((tmp[0],''))
	############--modified saling--############

	aln = Alignment(env)
	for (code, chain) in (code_and_chain):
	    mdl = Model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain))
	    aln.append_model(mdl, atom_files=code, align_codes=code+chain)

	for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
		                            ((1., 0.5, 1., 1., 1., 0.), False, True),
		                            ((1., 1., 1., 1., 1., 0.), True, False)):
	    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
		       rr_file='$(LIB)/as1.sim.mat', overhang=30,
		       gap_penalties_1d=(-450, -50),
		       gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
		       dendrogram_file='structureAlignment.tree',
		       alignment_type='tree', # If 'progresive', the tree is not
		                              # computed and all structues will be
		                              # aligned sequentially to the first
		       feature_weights=weights, # For a multiple sequence alignment only
		                                # the first feature needs to be non-zero
		       improve_alignment=True, fit=True, write_fit=write_fit,
		       write_whole_pdb=whole, output='ALIGNMENT QUALITY')

	aln.write(file='structureAlignment.pap', alignment_format='PAP')
	aln.write(file='structureAlignment.ali', alignment_format='PIR')

	aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
		   rr_file='$(LIB)/as1.sim.mat', overhang=30,
		   gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
		   gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',
		   alignment_type='progressive', feature_weights=[0]*6,
		   improve_alignment=False, fit=False, write_fit=True,
		   write_whole_pdb=False, output='QUALITY')

if __name__ == '__main__':
	main(args)
