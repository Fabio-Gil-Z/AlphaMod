import os
import glob
import json

def gdt_ts_calculator(path_to_gdt_ts_file, group_type, sequence, config):
	model_name = path_to_gdt_ts_file.split(sep='/')[-1]
	with open(f'{path_to_gdt_ts_file}', 'r') as file:
		data = file.read().splitlines()	

	for line, line_info in enumerate(data):
		if "NUMBER_OF_ATOMS_AA:" in line_info:
			CA_line_list = line_info.split()
			CA = CA_line_list[1]
		if "SUMMARY(GDT)" in line_info:
			GDT_line_list = line_info.split()
			GDT = GDT_line_list[6]

	for target, domain in config.items():
		sequenceConfig = target
		nt = domain
		if sequenceConfig in sequence:			
			try:
				# print(f'CA: {CA}')
				# print(f'GDT: {GDT}')
				# print(f'nt: {nt}')
				GDT_Score = (float(CA)*float(GDT)) / float(nt)
				GDT_Score = round(GDT_Score, 2)
				print(f"{target} - {os.path.basename(path_to_gdt_ts_file)[7:-4]} --> GDT Score: {GDT_Score}")
				print("\n")

				model_name_clean = model_name.replace("GDT_TS_", "")

				with open(f'{group_type}/results', 'a') as results_output:
					results_output.write(f'{sequence},{os.path.basename(group_type)},{model_name_clean[:-4]},{GDT_Score:.2f}\n')
			except:
				print(f'Error in file {path_to_gdt_ts_file}')
			break

def main():
	print("######################################################")
	print("#########--- CALCULATING GDT_TS SCORE ---#############")
	print("######################################################")
	cwd = os.getcwd()
	GDT_TS_output_folder = f'{cwd}/Structure_Assessment_Section/GDT_TS/GDT_TS_output_folder'

	get_list_of_sequences = glob.glob(f'{GDT_TS_output_folder}/*')
	get_list_of_sequences.sort()

	path_to_json = f'{cwd}/config.json'
	with open(f'{path_to_json}', 'r') as file:
		data = file.read()
		configuration = json.loads(data)
	config = configuration["AlphaFold_Prediction_List"]
	'''
	group_type, means each Target, i.e T1024-D1, T1025-D1, ... T1099-D1
	Will have a group type: 
		AlphaFold
		five_ranked_unsupervised
		2best_usupervised
		2best_supervised
	Inside each group_type there should be the GDT_TS files we need to calculate the GDT_TS score
	For example:
	GDT_TS_output_folder:
		T1024-D1:                       <-- Target
			-AlphaFold                  <-- group_type
				-GDT_TS_ranked_0.pdb
				-GDT_TS_ranked_1.pdb
				-GDT_TS_ranked_2.pdb
				-GDT_TS_ranked_3.pdb
				-GDT_TS_ranked_4.pdb
				-superposition_ranked_0.pdb
				-superposition_ranked_1.pdb
				-superposition_ranked_2.pdb
				-superposition_ranked_3.pdb
				-superposition_ranked_4.pdb
			-five_ranked_unsupervised   <-- group_type
			-2best_usupervised          <-- group_type
			-2best_supervised           <-- group_type
		T1025-D1:                       <-- Target
			-AlphaFold                  <-- group_type
			-five_ranked_unsupervised   <-- group_type
			-2best_usupervised          <-- group_type
			-2best_supervised           <-- group_type
	'''
	for sequence in get_list_of_sequences:
		get_group_type = glob.glob(f'{sequence}/*')
		if os.path.exists(f'{sequence}/results_all_types.csv'):
			os.remove(f'{sequence}/results_all_types.csv')
		for group_type in get_group_type:
			get_gdt_ts_files = glob.glob(f'{group_type}/GDT_TS_*.pdb')
			get_gdt_ts_files.sort()
			sequence_path = os.path.dirname(group_type)
			sequence_name = os.path.basename(sequence_path)
			if os.path.exists(f'{group_type}/results'):
				os.remove(f'{group_type}/results')
			for gdt_file in get_gdt_ts_files:
				gdt_ts_calculator(gdt_file, group_type, sequence_name, config)

if __name__=="__main__":
	main()