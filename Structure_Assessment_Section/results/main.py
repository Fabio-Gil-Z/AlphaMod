import os
import glob
import json
import pandas as pd

def readConfigJSON(cwd):
    path_to_json = f'{cwd}/config.json'
    with open(f'{path_to_json}', 'r') as file:
        data = file.read()
        configuration = json.loads(data)
    return configuration

def getResults(path_to_results_folder, measure_type, structure_assessment_path):
    Targets = glob.glob(f'{path_to_results_folder}/*')
    Targets.sort()
    final_results_path = f'{path_to_results_folder}/results.csv'
    if os.path.exists(final_results_path):
        os.remove(final_results_path)
    for target in Targets:
        Model_Types = glob.glob(f'{target}/*')
        Model_Types.sort()
        for model in Model_Types:
            with open(final_results_path, 'a') as output_file:
                if measure_type == "GDT_TS":
                    results_file_path = f'{path_to_results_folder}/{os.path.basename(target)}/{os.path.basename(model)}/results'
                    with open(results_file_path, 'r') as results_file:
                        data = results_file.read()
                    output_file.write(data)
                if measure_type == "MOLPROBITY":
                    results_file_path = f'{path_to_results_folder}/{os.path.basename(target)}/{os.path.basename(model)}/MOLPROBITY_results_raw_csv.csv'
                    with open(results_file_path, 'r') as results_file:
                        data = results_file.readlines()
                        for line in data[1:]:
                            model_name = line.split(sep=",")[0]
                            molprob_score = line.split(sep=",")[1]
                            output_file.write(f'{os.path.basename(target)},{os.path.basename(model)},{model_name[:-4]},{molprob_score}')
                if measure_type == "PROCHECK":
                    results_file_path = f'{path_to_results_folder}/{os.path.basename(target)}/{os.path.basename(model)}/PROCHECK_results_raw_csv.csv'
                    with open(results_file_path, 'r') as results_file:
                        data = results_file.readlines()
                        for line in data[1:]:
                            model_name = line.split(sep=",")[0]
                            procheck_score = line.split(sep=",")[1]
                            output_file.write(
                                f'{os.path.basename(target)},{os.path.basename(model)},{model_name[:-4]},{procheck_score}')
                if measure_type == "PROSA":
                    results_file_path = f'{path_to_results_folder}/{os.path.basename(target)}/{os.path.basename(model)}/PROSA_results_raw_csv.csv'
                    with open(results_file_path, 'r') as results_file:
                        data = results_file.readlines()
                        for line in data[1:]:
                            model_name = line.split(sep=",")[0]
                            PROSA_score = line.split(sep=",")[1]
                            output_file.write(
                                f'{os.path.basename(target)},{os.path.basename(model)},{model_name[:-4]},{PROSA_score}')
                if measure_type == "QMEAN":
                    results_file_path = f'{path_to_results_folder}/{os.path.basename(target)}/{os.path.basename(model)}/QMEAN_results_cleaned_csv.csv'
                    with open(results_file_path, 'r') as results_file:
                        data = results_file.readlines()
                        for line in data[1:]:
                            model_name = line.split(sep=",")[0]
                            QMEAN_score = line.split(sep=",")[1]
                            output_file.write(f'{os.path.basename(target)},{os.path.basename(model)},{model_name[:-4]},{QMEAN_score}')
    if measure_type == "GDT_TS":
        with open(f'{structure_assessment_path}/GDT_TS/GDT_TS_output_folder/results.csv', 'r') as file:
            data = file.readlines()
        with open(f'{structure_assessment_path}/GDT_TS/GDT_TS_output_folder/results.csv', 'w') as file:
            file.write("Target,Model Type,Model,GDT_TS\n")
            file.writelines(data[:])
    if measure_type == "QMEAN":
        with open(f'{structure_assessment_path}/QMEAN/results/results.csv', 'r') as file:
            data = file.readlines()
        with open(f'{structure_assessment_path}/QMEAN/results/results.csv', 'w') as file:
            file.write("Target,Model Type,Model,QMEAN\n")
            file.writelines(data[:])
    if measure_type == "PROSA":
        with open(f'{structure_assessment_path}/PROSA/results/results.csv', 'r') as file:
            data = file.readlines()
        with open(f'{structure_assessment_path}/PROSA/results/results.csv', 'w') as file:
            file.write("Target,Model Type,Model,PROSA\n")
            file.writelines(data[:])
    if measure_type == "MOLPROBITY":
        with open(f'{structure_assessment_path}/MOLPROBITY/results/results.csv', 'r') as file:
            data = file.readlines()
        with open(f'{structure_assessment_path}/MOLPROBITY/results/results.csv', 'w') as file:
            file.write("Target,Model Type,Model,MOLPROBITY\n")
            file.writelines(data[:])
    if measure_type == "PROCHECK":
        with open(f'{structure_assessment_path}/PROCHECK/results/results.csv', 'r') as file:
            data = file.readlines()
        with open(f'{structure_assessment_path}/PROCHECK/results/results.csv', 'w') as file:
            file.write("Target,Model Type,Model,PROCHECK\n")
            file.writelines(data[:])
    return pd.read_csv(final_results_path)

def bordaCalc(target_dataframe):
    df_plddt = target_dataframe.sort_values(by="pLDDT", ascending=False).reset_index(drop=True)
    df_Qmean = target_dataframe.sort_values(by="QMEAN", ascending=False).reset_index(drop=True)

    borda_score = {}
    for i in range(0, 5):
        borda_score_QMEAN = int(df_Qmean.loc[df_Qmean['Model'] == f"ranked_{i}"].index[0]) + 1
        borda_score_pLDDT = int(df_plddt.loc[df_plddt['Model'] == f"ranked_{i}"].index[0]) + 1
        borda_score_final = borda_score_QMEAN + borda_score_pLDDT
        borda_score[f'ranked_{i}.pdb'] = borda_score_final
    # print(borda_score)
    borda_score_list = [value for key, value in borda_score.items()]
    target_dataframe["BORDASCORE"] = borda_score_list
    return target_dataframe

def bordaScore(full_dataframe):
    scores_by_target = full_dataframe.groupby(full_dataframe.Target)
    tmp_dataframe = pd.DataFrame()
    for score in scores_by_target:
        alphafold_rows = score[1][score[1]["Model Type"] == "AlphaFold"].copy()
        new_score = bordaCalc(alphafold_rows)
        tmp_dataframe = pd.concat([tmp_dataframe, new_score])
    list_of_columns = full_dataframe.columns.tolist()
    full_dataframe = pd.merge(tmp_dataframe, full_dataframe, on=list_of_columns, how="right")
    return full_dataframe
def mergeDataFrame(df1, df2):
    final_df = pd.merge(df1, df2, on=["Target", "Model Type", "Model"], how="right")
    final_df = final_df.sort_values(by=["Target", "Model Type"])
    return final_df
def main():
    cwd = os.getcwd()
    structure_assessment_path = f'{cwd}/Structure_Assessment_Section'
    configuration = readConfigJSON(cwd)

    df_pLDDT = pd.read_csv(f'{structure_assessment_path}/pLDDT/pLDDT.csv')
    final_df = df_pLDDT
    df_QMEAN= getResults(f'{structure_assessment_path}/QMEAN/results/', "QMEAN",  structure_assessment_path)
    final_df = mergeDataFrame(df_QMEAN, df_pLDDT)
    final_df = bordaScore(final_df)

    order = ["Target", "Model Type", "Model", "GDT_TS", "pLDDT", "QMEAN", "BORDASCORE", "PROSA", "MOLPROBITY", "PROCHECK", "DOPESCORE"]
    column_list = final_df.columns.tolist()
    csv_column_save_order = []
    for column_name in order:
        if column_name in column_list:
            csv_column_save_order.append(column_name)

    if not configuration["StructureAssessment"]["first_run_flag"].upper() == "TRUE":
        final_df = final_df[csv_column_save_order]
        final_df.to_csv("results_summary.csv", index=False)
        os.remove(f'{structure_assessment_path}/results/results.csv')
    else:
        final_df.to_csv(f'{structure_assessment_path}/results/results.csv', index=False)

    '''
   Updating Modeller_Scripts/config.json file 
   Setting "first_run_flag" to ---> FALSE
   '''
    path_to_json = f'{cwd}/config.json'
    with open(f'{path_to_json}', 'r') as file:
        data = file.read()
    configuration = json.loads(data)
    configuration["StructureAssessment"]["first_run_flag"] = "FALSE"
    with open(f'{path_to_json}', 'w') as outputfile:
        json.dump(configuration, outputfile, indent=4)

if __name__=="__main__":
    main()