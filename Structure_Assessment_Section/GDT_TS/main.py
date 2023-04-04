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

def cleaning_MOL2_directory(GDT_TS_path):
    #Cleaning up MOL2 directory...
    MOL2_directory_path = f'{GDT_TS_path}/LGA_Zemla/MOL2/'
    to_clean_list = os.listdir(MOL2_directory_path)
    for pdb_to_delete in to_clean_list:
        if pdb_to_delete.endswith(".pdb"):
            os.remove(os.path.join(MOL2_directory_path, pdb_to_delete))
    #Cleaning up TMP directory...
    TMP_directory_path = f'{GDT_TS_path}/LGA_Zemla/TMP/'
    to_clean_list = os.listdir(TMP_directory_path)
    for file in to_clean_list:
        if not os.path.basename(file) == ".gitkeep":
            os.remove(os.path.join(TMP_directory_path, file))

def superposition_and_GDT(cwd, GDT_TS_path, Target_List, alphafold_or_modeller_results_path, output_folder):
    #################################################################
    ###                                                           ###
    ###     Script to get the superposition and GDT_TS files      ###
    ###                                                           ###
    #################################################################

    template_pdb_file_location = f'{cwd}/AlphaFold_Section/AlphaFold_Files/templates'
    template_pdbs_folders = glob.glob(f'{template_pdb_file_location}/*')
    template_pdbs_folders.sort()    

    '''
    We only calculate those Targets that are in the AlphaMod/config.json file --> AlphaFold_Prediction_List
    Target_List --> are the values in 'AlphaFold_Prediction_List'
    working_dict is a dictionary that has as key the Target and as values the paths to the pdbs
    For example, for Target T1024-D0 and T1024-D1 for ranked pdbs
    The working_dict would look like:  
        {"AlphaMod/AlphaFold_Section_AlphaFold_Files/templates/T1024-D1/T1024-D1.pdb": 
                  [AlphaMod/AlphaFold_Section_AlphaFold_Files/results/T1024-D0/ranked_0.pdb,
                   AlphaMod/AlphaFold_Section_AlphaFold_Files/results/T1024-D0/ranked_1.pdb,
                   AlphaMod/AlphaFold_Section_AlphaFold_Files/results/T1024-D0/ranked_2.pdb,
                   AlphaMod/AlphaFold_Section_AlphaFold_Files/results/T1024-D0/ranked_3.pdb,
                   AlphaMod/AlphaFold_Section_AlphaFold_Files/results/T1024-D0/ranked_4.pdb]
    where -> "AlphaMod/AlphaFold_Section_AlphaFold_Files/templates/T1024-D1/T1024-D1.pdb" is the KEY of the dictionary
    So we can use the 'superposition_and_GDT_TS.sh' script easier, because we need the ranked pdbs, but we always
    make the superposition to the same template.
    '''

    working_dict = {}
    for template in template_pdbs_folders:
        for key, value in Target_List.items():
            if key in template:
                for target in alphafold_or_modeller_results_path:
                    if key in target:
                        if output_folder == "AlphaFold":
                            working_dict[f'{template}/{os.path.basename(template)}.pdb'] = [pdb for pdb in glob.glob(f'{target}/ranked*.pdb')]
                        if output_folder == "two_best_supervised":
                            working_dict[f'{template}/{os.path.basename(template)}.pdb'] = [pdb for pdb in glob.glob(f'{target}/two_best_supervised/*.pdb')]
                        if output_folder == "two_best_unsupervised":
                            working_dict[f'{template}/{os.path.basename(template)}.pdb'] = [pdb for pdb in glob.glob(f'{target}/two_best_unsupervised/*.pdb')]
                        if output_folder == "five_ranked_unsupervised":
                            working_dict[f'{template}/{os.path.basename(template)}.pdb'] = [pdb for pdb in glob.glob(f'{target}/five_ranked_unsupervised/*.pdb')]
                        break

    # Creating Modeller_results folder
    GDT_TS_output_folder = f'{GDT_TS_path}/GDT_TS_output_folder'
    if not os.path.exists(GDT_TS_output_folder):
        os.mkdir(GDT_TS_output_folder)
        print("Directory: ", f'{GDT_TS_output_folder}', " -> has been Created ")

    for template_pdb_file_location, predicted_PDB_FILES_location in working_dict.items():
        predicted_PDB_FILES_location.sort()
        file_name_lenght = len(os.path.basename(template_pdb_file_location[:-4]))
        print("#"*56 + "#"*file_name_lenght)
        print("#"*3 + " "*5 + f"Calculating Superposition and GDT_TS of {os.path.basename(template_pdb_file_location[:-4])}" + " "*5 + "#"*3)
        print("#"*56 + "#"*file_name_lenght)
        # Creating folder for each template (e.g. T1024, T1024-D1, T1024-D2 ... T1099-D1)
        template_output_folder = f'{GDT_TS_output_folder}/{os.path.basename(template_pdb_file_location[:-4])}'
        if not os.path.exists(template_output_folder):
            os.mkdir(template_output_folder)
            print("Directory: ", f'{template_output_folder}', " -> has been Created ")
        # Creating folder for each AlphaFold
        gdt_ts_alphafold_output_folder = f'{template_output_folder}/{output_folder}'
        if not os.path.exists(gdt_ts_alphafold_output_folder):
            os.mkdir(gdt_ts_alphafold_output_folder)
            print("Directory: ", f'{gdt_ts_alphafold_output_folder}', " -> has been Created ")
        for PDB_FILE in predicted_PDB_FILES_location:
            subprocess.check_call(f"./superposition_and_GDT_TS.sh %s %s" % (PDB_FILE, template_pdb_file_location), shell=True, cwd=GDT_TS_path)
            print(f'    superposition_and_GDT_TS of {os.path.basename(PDB_FILE)} finished')
        print(f'\n')

        superposition_files = glob.glob(f'{GDT_TS_path}/LGA_Zemla/MOL2/GDT_TS_*.pdb')
        GDT_TS_files = glob.glob(f'{GDT_TS_path}/LGA_Zemla/MOL2/superposition_*.pdb')
        for superposition_file in superposition_files:
            shutil.copy(superposition_file, gdt_ts_alphafold_output_folder)
        for GDT_TS_file in GDT_TS_files:
            shutil.copy(GDT_TS_file, gdt_ts_alphafold_output_folder)
        # Cleaning up MOL2 directory...
        cleaning_MOL2_directory(GDT_TS_path)

def main():
    cwd = os.getcwd()
    configuration = readConfigJSON(cwd)
    Target_List = configuration["AlphaFold_Prediction_List"]
    first_run_flag = configuration["StructureAssessment"]["first_run_flag"].upper()
    GDT_TS_path = f'{cwd}/Structure_Assessment_Section/GDT_TS'   

    '''
    If first_run_flag == "TRUE", means we are running Structure_Assessment tools for the first time, which means
    we are going to calculate alphafold models only.
    '''
    if first_run_flag == "TRUE":
        alphafold_results_path = f'{cwd}/AlphaFold_Section/AlphaFold_Files/results'
        alphafold_Targets = glob.glob(f'{alphafold_results_path}/*')
        alphafold_Targets.sort()
        superposition_and_GDT(cwd, GDT_TS_path, Target_List, alphafold_Targets, "AlphaFold")
    else:
        run_two_best_supervised = configuration["StructureAssessment"]["2best_supervised"]
        run_two_best_unsupervised = configuration["StructureAssessment"]["2best_unsupervised"]
        run_five_ranked_unsupervised = configuration["StructureAssessment"]["5ranked_unsupervised"]
        modeller_results_path = f'{cwd}/Modeller_Section/results'
        modeller_Targets = glob.glob(f'{modeller_results_path}/*')
        modeller_Targets.sort()
        if run_two_best_supervised.upper() == "TRUE":
            superposition_and_GDT(cwd, GDT_TS_path, Target_List, modeller_Targets, "two_best_supervised")
        if run_two_best_unsupervised.upper() == "TRUE":
            superposition_and_GDT(cwd, GDT_TS_path, Target_List, modeller_Targets, "two_best_unsupervised")
        if run_five_ranked_unsupervised.upper() == "TRUE":
            superposition_and_GDT(cwd, GDT_TS_path, Target_List, modeller_Targets, "five_ranked_unsupervised")

if __name__=="__main__":
    main()