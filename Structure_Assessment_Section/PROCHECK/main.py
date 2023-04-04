import os
import glob
import json
import shutil
import subprocess

def readConfigJSON(cwd):
    path_to_json = f'{cwd}/config.json'
    with open(f'{path_to_json}', 'r') as file:
        data = file.read()
        configuration = json.loads(data)
    return configuration


def iterator(template_pdbs_folders, Target_List, alphafold_or_modeller_results_path, mode, output_folder, cwd, Google_Chrome_Selenium_Web_Bot_Folder, configuration):
    for template in template_pdbs_folders:
        for key, value in Target_List.items():
            if key in template:
                for target in alphafold_or_modeller_results_path:
                    if key in target:
                        if mode == "AlphaFold":
                            pdbs = glob.glob(f'{target}/ranked*.pdb')
                            pdbs.sort()
                        if mode == "two_best_supervised":
                            pdbs = glob.glob(f'{target}/{mode}/*.pdb')
                            pdbs.sort()
                        if mode == "two_best_unsupervised":
                            pdbs = glob.glob(f'{target}/{mode}/*.pdb')
                            pdbs.sort()
                        if mode == "five_ranked_unsupervised":
                            pdbs = glob.glob(f'{target}/{mode}/*.pdb')
                            pdbs.sort()
                        for pdb in pdbs:
                            shutil.copy(pdb, f'{Google_Chrome_Selenium_Web_Bot_Folder}/pdb_files')
                        PROCHECK = configuration["StructureAssessment"]["PROCHECK"]
                        while True:
                            subprocess.check_call(f'python3 {Google_Chrome_Selenium_Web_Bot_Folder}/main.py -Q false -P false -M false -C {PROCHECK} -T {target}', shell=True, cwd=cwd)
                            PROCHECK_result_files = glob.glob(f"{Google_Chrome_Selenium_Web_Bot_Folder}/PROCHECK/*")
                            # If everything went good, we will have 6 files and no errors otherwise we keep trying.
                            if len(PROCHECK_result_files) == 6:
                                break
                        if not os.path.exists(f'{output_folder}/{os.path.basename(target)}'):
                            os.mkdir(f'{output_folder}/{os.path.basename(target)}')
                        if not os.path.exists(f'{output_folder}/{os.path.basename(target)}/{mode}'):
                            os.mkdir(f'{output_folder}/{os.path.basename(target)}/{mode}')
                        
                        for file in PROCHECK_result_files:
                            shutil.move(file, f'{output_folder}/{os.path.basename(target)}/{mode}')
                        deletePDB()

def deletePDB():
    cwd = os.getcwd()
    pdb_files_folder = f'{cwd}/Structure_Assessment_Section/Google_Chrome_Selenium_Web_Bot/pdb_files'
    pdb_files_to_delete = glob.glob(f'{pdb_files_folder}/*')
    for pdb in pdb_files_to_delete:
        os.remove(pdb)
def PROCHECK(cwd, configuration):
    cwd = os.getcwd()
    Google_Chrome_Selenium_Web_Bot_Folder = f'{cwd}/Structure_Assessment_Section/Google_Chrome_Selenium_Web_Bot'
    configuration = readConfigJSON(cwd)
    template_pdb_file_location = f'{cwd}/AlphaFold_Section/AlphaFold_Files/templates'
    template_pdbs_folders = glob.glob(f'{template_pdb_file_location}/*')
    template_pdbs_folders.sort()
    alphafold_results_path = f'{cwd}/AlphaFold_Section/AlphaFold_Files/results'
    alphafold_Targets = glob.glob(f'{alphafold_results_path}/*')
    alphafold_Targets.sort()
    run_alphafold = configuration["StructureAssessment"]["AlphaFold_Models"]
    run_two_best_supervised = configuration["StructureAssessment"]["2best_supervised"]
    run_two_best_unsupervised = configuration["StructureAssessment"]["2best_unsupervised"]
    run_five_ranked_unsupervised = configuration["StructureAssessment"]["5ranked_unsupervised"]
    modeller_results_path = f'{cwd}/Modeller_Section/results'
    modeller_Targets = glob.glob(f'{modeller_results_path}/*')
    modeller_Targets.sort()
    configuration = readConfigJSON(cwd)
    Target_List = configuration["AlphaFold_Prediction_List"]
    output_folder = f'{cwd}/Structure_Assessment_Section/PROCHECK/results'
    if os.path.exists(f'{output_folder}/results.csv'):
        os.remove(f'{output_folder}/results.csv')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
        print("Directory: ", f'{output_folder}', " -> has been Created ")

    if run_alphafold.upper() == "TRUE":
        iterator(template_pdbs_folders, Target_List, alphafold_Targets, "AlphaFold", output_folder, cwd, Google_Chrome_Selenium_Web_Bot_Folder, configuration)
    if run_two_best_supervised.upper() == "TRUE":
        iterator(template_pdbs_folders, Target_List, modeller_Targets, "two_best_supervised", output_folder, cwd, Google_Chrome_Selenium_Web_Bot_Folder, configuration)
    if run_two_best_unsupervised.upper() == "TRUE":
        iterator(template_pdbs_folders, Target_List, modeller_Targets, "two_best_unsupervised", output_folder, cwd, Google_Chrome_Selenium_Web_Bot_Folder, configuration)
    if run_five_ranked_unsupervised.upper() == "TRUE":
        iterator(template_pdbs_folders, Target_List, modeller_Targets, "five_ranked_unsupervised", output_folder, cwd, Google_Chrome_Selenium_Web_Bot_Folder, configuration)
    os.rmdir(f"{Google_Chrome_Selenium_Web_Bot_Folder}/PROCHECK")

def main():
    cwd = os.getcwd()
    configuration = readConfigJSON(cwd)
    PROCHECK(cwd, configuration)

if __name__=="__main__":
    main()