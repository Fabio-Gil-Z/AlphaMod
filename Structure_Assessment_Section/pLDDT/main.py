import os
import glob
import json

def readPLDDT(path):
    with open(f'{path}', 'r') as file:
        data = file.read()
    return json.loads(data)
def readConfigJSON(cwd):
    path_to_json = f'{cwd}/config.json'
    with open(f'{path_to_json}', 'r') as file:
        data = file.read()
        configuration = json.loads(data)
    return configuration

def main():
    cwd = os.getcwd()

    fasta_files_path = f'{cwd}/AlphaFold_Section/AlphaFold_Files/Fasta_Files/'
    fasta_files = glob.glob(f'{fasta_files_path}/*')
    fasta_files.sort()


    alphafold_results_path = f'{cwd}/AlphaFold_Section/AlphaFold_Files/results'
    alphafold_Targets = glob.glob(f'{alphafold_results_path}/*')
    alphafold_Targets.sort()

    configuration = readConfigJSON(cwd)
    Target_List = configuration["AlphaFold_Prediction_List"]
    output_folder = f'{cwd}/Structure_Assessment_Section/pLDDT/'
    if os.path.exists(f'{output_folder}/pLDDT.csv'):
        os.remove(f'{output_folder}/pLDDT.csv')


    for fasta_file in fasta_files:
        for key, value in Target_List.items():
            if key in fasta_file:
                for target in alphafold_Targets:
                    if key in target:
                        plddt_file = f'{target}/ranking_debug.json'
                        plddt_data = readPLDDT(plddt_file)
                        tmp_dict = {}
                        ranked_counter = 0
                        with open(f'{output_folder}/pLDDT.csv', 'a') as output_file:
                            for ranked in plddt_data["order"]:
                                tmp_dict[f'ranked_{ranked_counter}'] = plddt_data["plddts"][ranked]
                                ranked_counter += 1
                            for model, plddt in tmp_dict.items():
                                output_file.write(f"{os.path.basename(target)},AlphaFold,{model},{round(float(plddt), 2)}\n")

    with open(f'{output_folder}/pLDDT.csv', 'r') as file:
        data = file.readlines()
    with open(f'{output_folder}/pLDDT.csv', 'w') as output_file:
        output_file.write("Target,Model Type,Model,pLDDT\n")
        output_file.writelines(data[:])
if __name__ == "__main__":
    main()