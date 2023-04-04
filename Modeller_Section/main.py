import os
import glob
import json
import shutil
import subprocess
import pandas as pd

def readConfigJSON(cwd):
    path_to_json = f'{cwd}/config.json'
    with open(f'{path_to_json}', 'r') as file:
        data = file.read()
        configuration = json.loads(data)
    return configuration
def cleaner(directory_to_clean):
    list_of_files_we_do_not_want_to_delete = ['align2d_mult.py', 'fastaToAli.py', 'model_mult.py', 'salign.py']
    list_of_all_files = glob.glob(f'{directory_to_clean}/*')
    for file_to_delete in list_of_all_files:
        file_to_delete_basename = os.path.basename(file_to_delete)
        if not file_to_delete_basename in list_of_files_we_do_not_want_to_delete:
            os.remove(file_to_delete)
def runModellerScripts(cwd, output_folder, mode):
    modeller_scripts_folder_path = f'{cwd}/Modeller_Section/Modeller_Scripts'
    modeller_without_installation_path = f'{cwd}/Modeller_Section/modeller10.4/bin/modpy.sh'
    fasta_file = glob.glob(f'{modeller_scripts_folder_path}/*.fasta')[0]
    print(fasta_file)
    models = glob.glob(f'{modeller_scripts_folder_path}/*.pdb')
    models_basename = []

    '''
    models_basename is necessary to send it as input to saling.py and model_mult.py
    '''

    for model in models:
        models_basename.append(f'{os.path.basename(model)}')
    subprocess.check_call(f'python3 fastaToAli.py -d {fasta_file}', shell=True, cwd=modeller_scripts_folder_path)
    subprocess.check_call(f'{modeller_without_installation_path} python3 salign.py -l {models_basename}', shell=True, cwd=modeller_scripts_folder_path)
    subprocess.check_call(f'{modeller_without_installation_path} python3 align2d_mult.py', shell=True, cwd=modeller_scripts_folder_path)
    subprocess.check_call(f'{modeller_without_installation_path} python3 model_mult.py -l {models_basename}', shell=True, cwd=modeller_scripts_folder_path)
    modeller_models_files_paths = glob.glob(f'{modeller_scripts_folder_path}/inputSequence*.pdb')
    if mode == "five_ranked":
        for i, modeller_model in enumerate(modeller_models_files_paths):
            shutil.copy(modeller_model, f'{output_folder}/Mod-5ranked_unsupervised_{i}.pdb')
    else:
        for i, modeller_model in enumerate(modeller_models_files_paths):
            shutil.copy(modeller_model, f'{output_folder}/Mod-2best_{mode}_{i}.pdb')
    '''
    cleaning modeller_scripts directory
    '''
    cleaner(modeller_scripts_folder_path)

def two_best(cwd, results_folder_path, mode):
    '''
    two_best_supervised is based on GDT_TS
    We will select the two_best models based on GDT_TS and then give them as input to Modeller.
    We want to select the models with the HIGHEST GDT_TS score

    two_best_unsupervised is based on pLDDT and QMEAN.
    Using the BORDA COUNT voting system we can obtain BORDASCORE.
    BORDASCORE will be used to select the best two models and then give them as input to Modeller.
    We want to select the models with the LOWEST BORDASCORE score
    '''
    if mode == "supervised":
        sort_mode = "GDT_TS"
        ascending_mode = False
    elif mode == "unsupervised":
        sort_mode = "BORDASCORE"
        ascending_mode = True
    else:
        print("Mode not recognized at: Modeller_Section/main.py\nFunction: two_best() line: 46.\nPlease select one of the following: ['supervised', 'unsupervised']")

    results_file_path = f'{cwd}/Structure_Assessment_Section/results/results.csv'
    alphafold_fasta_files_folder_path = f'{cwd}/AlphaFold_Section/AlphaFold_Files/Fasta_Files'
    alphafold_results_folder_path = f'{cwd}/AlphaFold_Section/AlphaFold_Files/results'
    modeller_scripts_folder_path = f'{cwd}/Modeller_Section/Modeller_Scripts'
    df = pd.read_csv(results_file_path)
    df_grouped_by_targets = df.groupby("Target")
    if os.path.exists(f'{results_folder_path}/List_of_2best_models_selected._based_on_{sort_mode}csv'):
        os.remove(f'{results_folder_path}/List_of_2best_models_selected_based_on_{sort_mode}.csv')
    with open(f'{results_folder_path}/List_of_2best_models_selected_based_on_{sort_mode}.csv', 'a') as outputfile:
        outputfile.write("Target,Model\n")
    for target_df in df_grouped_by_targets:
        sorted_target_df = target_df[1].sort_values(by=[sort_mode], ascending=ascending_mode, ignore_index=True)
        target = sorted_target_df.iloc[0]["Target"]
        Model_Type = f'two_best_{mode}'

        best_model_1 = sorted_target_df.iloc[0]
        best_model_2 = sorted_target_df.iloc[1]
        if not os.path.exists(f'{results_folder_path}/{target}'):
            os.mkdir(f'{results_folder_path}/{target}')
        if not os.path.exists(f'{results_folder_path}/{target}/{Model_Type}'):
            os.mkdir(f'{results_folder_path}/{target}/{Model_Type}')
        output_folder = f'{results_folder_path}/{target}/{Model_Type}'

        with open(f'{results_folder_path}/List_of_2best_models_selected_based_on_{sort_mode}.csv', 'a') as outputfile:
            outputfile.write(f'{best_model_1["Target"]},{best_model_1["Model"]}\n')
            outputfile.write(f'{best_model_2["Target"]},{best_model_2["Model"]}\n')

        # Copy corresponding target fasta file to --> modeller_scripts_folder_path
        shutil.copy(f'{alphafold_fasta_files_folder_path}/{target}.fasta', f'{modeller_scripts_folder_path}')

        # Copy the selected two best targets files to --> modeller_scripts_folder_path
        shutil.copy(f'{alphafold_results_folder_path}/{target}/{best_model_1["Model"]}.pdb', f'{modeller_scripts_folder_path}/mod1.pdb')
        shutil.copy(f'{alphafold_results_folder_path}/{target}/{best_model_2["Model"]}.pdb', f'{modeller_scripts_folder_path}/mod2.pdb')

        runModellerScripts(cwd, output_folder, mode)
def five_ranked_unsupervised(cwd, results_folder_path):
    results_file_path = f'{cwd}/Structure_Assessment_Section/results/results.csv'
    alphafold_fasta_files_folder_path = f'{cwd}/AlphaFold_Section/AlphaFold_Files/Fasta_Files'
    alphafold_results_folder_path = f'{cwd}/AlphaFold_Section/AlphaFold_Files/results'
    modeller_scripts_folder_path = f'{cwd}/Modeller_Section/Modeller_Scripts'
    df = pd.read_csv(results_file_path)
    df_grouped_by_targets = df.groupby("Target")
    for target_df in df_grouped_by_targets:
        target = target_df[1].iloc[0]["Target"]
        Model_Type = "five_ranked_unsupervised"

        if not os.path.exists(f'{results_folder_path}/{target}'):
            os.mkdir(f'{results_folder_path}/{target}')
        if not os.path.exists(f'{results_folder_path}/{target}/{Model_Type}'):
            os.mkdir(f'{results_folder_path}/{target}/{Model_Type}')
        output_folder = f'{results_folder_path}/{target}/{Model_Type}'

        # Copy corresponding target fasta file to --> modeller_scripts_folder_path
        shutil.copy(f'{alphafold_fasta_files_folder_path}/{target}.fasta', f'{modeller_scripts_folder_path}')

        # Copy 5 ranked models to --> modeller_scripts_folder_path
        for i in range(5):
            shutil.copy(f'{alphafold_results_folder_path}/{target}/ranked_{i}.pdb',
                        f'{modeller_scripts_folder_path}/mod{i}.pdb')

        runModellerScripts(cwd, output_folder, "five_ranked")

def main():
    cwd = os.getcwd()
    results_folder_path = f'{cwd}/Modeller_Section/results'
    configuration = readConfigJSON(cwd)
    two_best_supervised_flag = configuration["StructureAssessment"]["2best_supervised"]
    GDT_TS_flag = configuration["StructureAssessment"]["GDT_TS"]
    two_best_unsupervised_flag = configuration["StructureAssessment"]["2best_unsupervised"]
    five_ranked_unsupervised_flag = configuration["StructureAssessment"]["5ranked_unsupervised"]

    modlib_modeller_config_file = f'{cwd}/Modeller_Section/modeller10.4/modlib/modeller/config.py'
    install_dir = f"install_dir = r'{cwd}/Modeller_Section/modeller10.4/'"
    with open(modlib_modeller_config_file, 'r') as read_file:
        data = read_file.readlines()
    data[0] = install_dir
    with open(modlib_modeller_config_file, 'w') as write_file:
        for line in data:
            write_file.write(line + "\n")

    if two_best_supervised_flag.upper() == "TRUE" and GDT_TS_flag.upper() == "TRUE":
        two_best(cwd, f'{results_folder_path}', "supervised")
    if two_best_unsupervised_flag.upper() == "TRUE":
        two_best(cwd, f'{results_folder_path}', "unsupervised")
    if five_ranked_unsupervised_flag.upper() == "TRUE":
        five_ranked_unsupervised(cwd, f'{results_folder_path}')

if __name__ == '__main__':
    main()