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
    if configuration["Install_Requirements"] == "True".upper():
        print_config = json.dumps(configuration, indent=4)
        print("#" * 49)
        print("#" * 3 + " " * 15 + "CONFIGURATION" + " " * 15 + "#" * 3)
        print("#" * 49)
        print(print_config)
        print("First time launching AlphaMod, requirements will be installed...")
        input("Please press Enter...")
        subprocess.check_call(f'pip3 install -r {cwd}/requirements.txt', shell=True, cwd=cwd)
        configuration["Install_Requirements"] = "Already Installed"
        with open(f'{path_to_json}', 'w') as outputfile:
            json.dump(configuration, outputfile, indent=4)
        print("\n\n")
        print("#" * 121)
        print("#" * 3 + " " * 15 + "Requirements have been installed, press enter to continue, AlphaFold will be launched" + " " * 15 + "#" * 3)
        input("#" * 121)
        return configuration
    else:
        print("#" * 49)
        print("#" * 3 + " " * 15 + "CONFIGURATION" + " " * 15 + "#" * 3)
        print("#" * 49)
        print_config = json.dumps(configuration, indent=4)
        print(print_config)
        return configuration
    
def pdbCleaner(folderPath):
    # The reason behind cleaning the templates is to be able to use them in the 4 tools automation
    # These are 4 different websites with different requirements
    # If we do not leave the pdbs with only the ATOM part, we will find that the website does not
    # Make the calculations and answers with an error on the template.
    # The websites:
    # QMEAN: https://swissmodel.expasy.org/qmean/
    # ProsaWeb: https://prosa.services.came.sbg.ac.at/prosa.php
    # MolProbity: http://molprobity.biochem.duke.edu/
    # Procheck: https://saves.mbi.ucla.edu/

    templateList = glob.glob(f'{folderPath}/*.pdb')
    templateList.sort()

    for template in templateList:
        new_pdb = []
        with open(template, 'r') as data:
            pdb = data.readlines()

        for line in pdb:
            line_as_list = line.split()
            # print(line_as_list, len(line_as_list))
            try:
                if line_as_list[0][0] == "A" and \
                   line_as_list[0][1] == "T" and \
                   line_as_list[0][2] == "O" and \
                   line_as_list[0][3] == "M":
                       new_pdb.append(line)
            except:
                pass

        with open(template, 'w') as data:
            data.writelines(new_pdb[:])
            
def runAlphaFold_Section(cwd):
    multiple_runs_main_path = f'{cwd}/AlphaFold_Section/AlphaFold_Multiple_Runs/main.py'
    subprocess.check_call(f'python3 {multiple_runs_main_path}', shell=True, cwd=cwd)
    if not os.path.exists(f'{cwd}/AlphaFold_Section/logs'):
        os.mkdir(f'{cwd}/AlphaFold_Section/logs')
        print(f"Folder 'logs' has been created at {cwd}/AlphaFold_Section/logs")
        print("All AlphaFold prediction logs have been stored in 'logs' folder.")
    log_list = glob.glob(f'{cwd}/*.log')
    for log in log_list:
        shutil.move(log, f'{cwd}/AlphaFold_Section/logs')

def runStructureAssessment(cwd, configuration):
    structure_assessment_path = f'{cwd}/Structure_Assessment_Section/'
    template_pdb_file_location = f'{cwd}/AlphaFold_Section/AlphaFold_Files/templates'
    template_pdbs_folders = glob.glob(f'{template_pdb_file_location}/*')
    template_pdbs_folders.sort()
    alphafold_results_path = f'{cwd}/AlphaFold_Section/AlphaFold_Files/results'
    alphafold_Targets = glob.glob(f'{alphafold_results_path}/*')
    alphafold_Targets.sort()
    pdbCleaner(template_pdbs_folders)
    pdbCleaner(alphafold_Targets)

    # pLDDT
    if configuration["StructureAssessment"]["first_run_flag"].upper() == "TRUE":
        if configuration["StructureAssessment"]["pLDDT"].upper() == "TRUE":
            subprocess.check_call(f'python3 {structure_assessment_path}/pLDDT/main.py', shell=True, cwd=cwd)
    
    # GDT_TS
    if configuration["StructureAssessment"]["GDT_TS"].upper() == "TRUE" and configuration["StructureAssessment"]["2best_supervised"].upper() == "TRUE":
        print("###################################################################################################")
        print("###                                                                                             ###")
        print("###                     To be able to calculate the GDT_TS Score                                ###")
        print("###                   we need to give permission to Zemla's scripts                             ###")
        print("###                         please enter your 'sudo' password                                   ###")
        print("###                        if no password is asked, is because                                  ###")
        print("###                           you already gave the permission                                   ###")
        print("###                           and this message can be ignored                                   ###")
        print("###                      The files we are giving permission are:                                ###")
        print("###       AlphaMod/Structure_Assessment_Section/GDT_TS/Zemla_GDT_TS_admin_permissions.sh        ###")
        print("###       AlphaMod/Structure_Assessment_Section/GDT_TS/superposition_and_GDT_TS.sh              ###")
        print("###       AlphaMod/Structure_Assessment_Section/GDT_TS/LGA_Zemla/lga                            ###")
        print("###       AlphaMod/Structure_Assessment_Section/GDT_TS/LGA_Zemla/MOL2/collect_PDB.pl            ###")
        print("###                                                                                             ###")
        print("###################################################################################################")
        GDT_TS_path = f'{cwd}/Structure_Assessment_Section/GDT_TS'
        subprocess.check_call("./Zemla_GDT_TS_admin_permissions.sh", shell=True, cwd=f'{GDT_TS_path}')
        subprocess.check_call(f'python3 {structure_assessment_path}/GDT_TS/main.py', shell=True, cwd=cwd)
        subprocess.check_call(f'python3 {structure_assessment_path}/GDT_TS/GDT_TS_Score_Calculator.py', shell=True, cwd=f'{cwd}')

    # QMEAN
    if configuration["StructureAssessment"]["QMEAN"].upper() == "TRUE":
        subprocess.check_call(f'python3 {structure_assessment_path}/QMEAN/main.py', shell=True, cwd=cwd)

    # PROSA
    if not configuration["StructureAssessment"]["first_run_flag"].upper() == "TRUE":
        if configuration["StructureAssessment"]["PROSA"].upper() == "TRUE":
            subprocess.check_call(f'python3 {structure_assessment_path}/PROSA/main.py', shell=True, cwd=cwd)

    # MOLPROBITY
    if not configuration["StructureAssessment"]["first_run_flag"].upper() == "TRUE":
        if configuration["StructureAssessment"]["MOLPROBITY"].upper() == "TRUE":
            subprocess.check_call(f'python3 {structure_assessment_path}/MOLPROBITY/main.py', shell=True, cwd=cwd)

    # PROCHECK
    if not configuration["StructureAssessment"]["first_run_flag"].upper() == "TRUE":
        if configuration["StructureAssessment"]["PROCHECK"].upper() == "TRUE":
            subprocess.check_call(f'python3 {structure_assessment_path}/PROCHECK/main.py', shell=True, cwd=cwd)

    # DOPESCORE
    if not configuration["StructureAssessment"]["first_run_flag"].upper() == "TRUE":
        if configuration["StructureAssessment"]["DOPESCORE"].upper() == "TRUE":
            subprocess.check_call(f'python3 {structure_assessment_path}/DOPESCORE/main.py', shell=True, cwd=cwd)

    # RESULTS
    subprocess.check_call(f'python3 {structure_assessment_path}/results/main.py', shell=True, cwd=cwd)

def runModeller_Section(cwd):
    modeller_main_file_path = f'{cwd}/Modeller_Section/main.py'
    subprocess.check_call(f'python3 {modeller_main_file_path}', shell=True, cwd=cwd)
def main():
    cwd = os.getcwd()
    try:
        os.environ["LD_LIBRARY_PATH"] = f"{os.environ['LD_LIBRARY_PATH']}:{cwd}/Modeller_Section/modeller10.4/lib/x86_64-intel8"
    except:
        os.environ["LD_LIBRARY_PATH"] = f"{cwd}/Modeller_Section/modeller10.4/lib/x86_64-intel8"
    try:
        os.environ["PYTHONPATH"] = f"{os.environ['PYTHONPATH']}:{cwd}/Modeller_Section/modeller10.4/modlib"
    except:
        os.environ["PYTHONPATH"] = f"{cwd}/Modeller_Section/modeller10.4/modlib"
    os.environ["PYTHONPATH"] = f"{os.environ['PYTHONPATH']}:{cwd}/Modeller_Section/modeller10.4/lib/x86_64-intel8/python3.3"


    configuration = readConfigJSON(cwd)
    # runAlphaFold_Section(cwd)
    runStructureAssessment(cwd, configuration)
    configuration = readConfigJSON(cwd)
    runModeller_Section(cwd)
    runStructureAssessment(cwd, configuration)


if __name__=="__main__":
    main()