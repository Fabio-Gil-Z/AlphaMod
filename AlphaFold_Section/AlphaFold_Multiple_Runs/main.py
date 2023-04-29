import os
import glob
import json
import subprocess

def readConfigJSON(cwd):
    path_to_json = f'{cwd}/config.json'
    with open(f'{path_to_json}', 'r') as file:
        data = file.read()
        configuration = json.loads(data)
    return configuration["AlphaFold_Databases"]
def readingTargetList():
	'''
	Reads the file "config.json" located at folder: AlphaMod/config.json
	Returns: a list with the names of the Targets written in the file "AlphaMod/config.json"
	'''
	print("############--- AlphaMod/AlphaFold_Section/AlphaFold_Multiple_Runs/main.py ---##############")
	print("#################--- function 'readingTargetList()' ---#####################")
	with open("config.json", 'r') as config_json_file:
		config_json_file_data = config_json_file.read()
		configuration = json.loads(config_json_file_data)
	list_of_targets = [target for target in configuration["AlphaFold_Prediction_List"]]
	print("Your Target List is:")
	for target in list_of_targets:
		print(target)
	print(f'The number of Targets entered as input are: {len(list_of_targets)}')
	print("##############################################################################")
	print("\n")
	list_of_targets.sort()
	return list_of_targets


def getFastaFilesPaths(cwd, list_of_targets):
	'''
	Receives as input "list_of_targets"
	Returns: a list with the paths to the fasta files to be calculated by AlphaFold
	'''
	print("############--- AlphaMod/AlphaFold_Section/AlphaFold_Multiple_Runs/main.py ---##############")
	print("#################--- function 'getFastaFilesPaths()' ---######################")
	paths_of_FastaFiles_to_Calculate = []

	path_to_alphafold_fasta_files = f'{cwd}/AlphaFold_Section/AlphaFold_Files/Fasta_Files'
	list_of_FastaFiles = glob.glob(f'{path_to_alphafold_fasta_files}/*')
	list_of_FastaFiles.sort()

	'''
	Independently of the number of fasta files stored at AlphaMod/AlphaFold_Section/AlphaFold_Files/Fasta_Files
	AlphaFold_multiple_runs/main.py will only run the Targets written in AlphaMod/config.json -> 'AlphaFold_Prediction_List' 
	'''
	for target in list_of_targets:
		for FastaFile in list_of_FastaFiles:
			if target in FastaFile:
				paths_of_FastaFiles_to_Calculate.append(FastaFile)
	paths_of_FastaFiles_to_Calculate.sort()
	print("The paths of the Fasta Files to calculate are:")
	for fastaFile_path in paths_of_FastaFiles_to_Calculate:
		print(fastaFile_path)
	print("##############################################################################")
	print("\n")
	return paths_of_FastaFiles_to_Calculate


def runAlphaFold(cwd, fastaFile_path, alphafold_params):
	fasta_file_basename = os.path.basename(fastaFile_path)
	print(f'#####################-- RUNNING {fasta_file_basename} --#########################################################')
	print(f'File location: {fastaFile_path}\n\n')
	fasta_file_basename = fasta_file_basename.split(sep='.')
	fasta_file_basename = fasta_file_basename[0]

	path_to_run_docker = f'{cwd}/AlphaFold_Section/alphafold/docker/run_docker.py'
	calling_template = f'python3 {path_to_run_docker} \
	--fasta_paths={alphafold_params["fasta_paths"]}/{fasta_file_basename}.fasta \
	--max_template_date={alphafold_params["max_template_date"]} \
	--data_dir={alphafold_params["data_dir"]} \
	--output_dir={alphafold_params["output_dir"]} \
	--logtostderr 2>&1 | tee {fasta_file_basename}.log'

	subprocess.check_call(f'{calling_template}', shell=True, cwd=f'{cwd}')



def launchingAlphaFold(cwd, paths_of_FastaFiles_to_Calculate):
	'''
	Receives as input "paths_of_FastaFiles_to_Calculate"

	AlphaFold calling template: found at 'https://github.com/deepmind/alphafold'
	python3 docker/run_docker.py \
	  --fasta_paths=your_protein.fasta \
	  --max_template_date=2022-01-01 \
	  --data_dir=$DOWNLOAD_DIR \
	  --output_dir=/home/user/absolute_path_to_the_output_dir
	'''

	print("############--- AlphaMod/AlphaFold_Section/AlphaFold_Multiple_Runs/main.py ---################")
	print("#################--- function 'launchingAlphaFold()' ---######################")


	default_fasta_paths = f'{cwd}/AlphaFold_Section/AlphaFold_Files/Fasta_Files'
	default_max_template_date = "2020-05-14"
	default_data_dir = readConfigJSON(cwd)
	default_output_dir = f'{cwd}/AlphaFold_Section/AlphaFold_Files/results'


	alphafold_params = {"fasta_paths":f'{default_fasta_paths}',
						 "max_template_date":f'{default_max_template_date}',
						 "data_dir":f'{default_data_dir}',
						 "output_dir":f'{default_output_dir}'}
		
	print("\n")
	print("################# AlphaFold Parameters #################")
	for key, value in alphafold_params.items():
		print(f'{key}: ', value)
	print("################# AlphaFold Parameters #################")
	print("\n")

	'''
	Running DeepMind's AlphaFold using as input the Fasta Files found at: AlphaMod/AlphaFold_Section/AlphaFold_Files/Fasta_Files
	Keep in mind, that it will only calculate those Targets given as input at: AlphaMod/config.json
	'''
	for fastaFile_path in paths_of_FastaFiles_to_Calculate:
		runAlphaFold(cwd, fastaFile_path, alphafold_params)
		print("\n\n")

def main():
	cwd = os.getcwd()
	list_of_targets = readingTargetList()
	paths_of_FastaFiles_to_Calculate = getFastaFilesPaths(cwd, list_of_targets)
	launchingAlphaFold(cwd, paths_of_FastaFiles_to_Calculate)

if __name__=="__main__":
	main()
