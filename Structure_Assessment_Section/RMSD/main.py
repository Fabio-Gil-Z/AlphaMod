import os
import glob
import subprocess

import pandas as pd


def writingRMSD(protein, tmp_log_path, RMSD_score_path):
    with open(tmp_log_path, 'r') as file:
        file_text = file.readlines()
    for i in file_text:
        print(i)
        if "pdb_file_name_moving =" in i:
            pdb_name = i.split()[2]
            print(pdb_name)
        if "RMSD between fixed and moving atoms (final):" in i:
            rmsd_score = i.split()[-1]
            print(rmsd_score)
            print(pdb_name)
            if "ranked" in pdb_name:
                if "Mod-5ranked_" in pdb_name:
                    with open(RMSD_score_path, 'a') as file:
                        file.write(f'{protein},Mod-5ranked,{pdb_name[:-4]},{rmsd_score}\n')
                else:
                    with open(RMSD_score_path, 'a') as file:
                        file.write(f'{protein},ranked,{pdb_name[:-4]},{rmsd_score}\n')
            if "2best" in pdb_name:
                if "Mod-2best_unsupervised_" in pdb_name:
                    with open(RMSD_score_path, 'a') as file:
                        file.write(f'{protein},Modalpha_2best_unsupervised,{pdb_name[:-4]},{rmsd_score}\n')
                else:
                    with open(RMSD_score_path, 'a') as file:
                        file.write(f'{protein},Modalpha_2best,{pdb_name[:-4]},{rmsd_score}\n')


def RMSD_calculator(folder_path):
    protein_chain_list = {"7MBY-D1": "B", "7ME0-D1": "A", "7EV9-D1": "A", "7LS5-D1": "A", "7EDA-D1": "A",
                          "7LCI-D1": "R", "7LVR-D1": "A", "7C2K-D1": "A", "7M7B-D1": "A", "7N8I-D1": "L",
                          "7MJS-D1": "H", "7L1K-D1": "A", "7L6U-D1": "A", "7KU7-D1": "A", "7KZZ-D1": "B",
                          "7LX5-D1": "B", "7BRM-D1": "A", "7LSX-D1": "A", "7LC6-D1": "A", "7MLZ-D1": "A",
                          "7MSW-D1": "A", "7RB9-D1": "B", "7BXT-D1": "A", "7M9C-D1": "A", "7LV9-D1": "B"}
    protein = os.path.basename(folder_path)
    if protein in protein_chain_list:
        print(f"Protein: {protein}, Chain: {protein_chain_list[protein]}")
        chain = protein_chain_list[protein]

    RMSD_score_path = f'{folder_path}/RMSD_scores.csv'
    protein_name = folder_path.split(sep="/")[-1]
    print(protein_name)
    if os.path.exists(RMSD_score_path):
        os.remove(RMSD_score_path)
    for i in range(5):
        subprocess.check_call(f"phenix.superpose_pdbs {protein_name}.pdb ranked_{i}.pdb selection_fixed=\"chain {chain} and name CA\" selection_moving=\"chain A and name CA\" > log_tmp.txt", shell=True, cwd=folder_path)
        writingRMSD(protein, f'{folder_path}/log_tmp.txt', RMSD_score_path)
        subprocess.check_call(f"phenix.superpose_pdbs {protein_name}.pdb Mod-5ranked_unsupervised_{i}.pdb selection_fixed=\"chain {chain} and name CA\" selection_moving=\"chain A and name CA\" > log_tmp.txt", shell=True, cwd=folder_path)
        writingRMSD(protein, f'{folder_path}/log_tmp.txt', RMSD_score_path)
        subprocess.check_call(f"phenix.superpose_pdbs {protein_name}.pdb Mod-2best_unsupervised_{i}.pdb selection_fixed=\"chain {chain} and name CA\" selection_moving=\"chain A and name CA\" > log_tmp.txt", shell=True, cwd=folder_path)
        writingRMSD(protein, f'{folder_path}/log_tmp.txt', RMSD_score_path)

    os.remove(f'{folder_path}/log_tmp.txt')
    # input(f"Ok Finished -> {protein_name}")

def final_csv(final_csv_path, Terwilliger_et_al_proteins):
    for target in Terwilliger_et_al_proteins:
        tmp_df = pd.read_csv(target)
        print(tmp_df)
        input()
def main():
    cwd = os.getcwd()
    Terwilliger_et_al_proteins = glob.glob(f'{cwd}/Structure_Assessment_Section/RMSD/phenix_superpose_pdbs/*')
    Terwilliger_et_al_proteins.sort()

    final_csv_path = f'{cwd}/summary_rmsd.csv'
    if os.path.exists(final_csv_path):
        os.remove(final_csv_path)

    for folder in Terwilliger_et_al_proteins:
        print(folder)
        try:
            RMSD_calculator(folder)
            with open(f"{folder}/RMSD_scores.csv", 'r') as file:
                info = file.readlines()
            with open(final_csv_path, 'a') as final_file:
                final_file.writelines(info[:])
        except:
            exit()

if __name__ == '__main__':
    main()


