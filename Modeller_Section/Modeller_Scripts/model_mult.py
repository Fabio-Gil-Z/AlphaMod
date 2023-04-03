from modeller import *
from modeller.automodel import *

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
	############--modified model_mult--############
	working_pdb_list = []
	for pdb in args.templateList:
		working_pdb_list.append(f'{cleaning_parser_list(pdb)}')
		
	code_and_chain = []
	for template in working_pdb_list:
    		code_and_chain.append(template.replace('_',''))
	############--modified model_mult--############
	env = Environ()
	a = AutoModel(env, alnfile='inputSequence-mult.ali',
		      knowns=(code_and_chain), sequence='inputSequence', assess_methods=(assess.DOPE,assess.GA341))
	a.starting_model = 1
	a.ending_model = 5
	a.make()
	
if __name__ == '__main__':
	main(args)
