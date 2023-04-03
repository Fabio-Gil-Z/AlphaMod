#!/bin/bash
################################ IMPORTANT ##############################
#####										                                            #####
#####  Do not delete TMP folder we need it, Zemla algorithm uses it #####
#####										                                            #####
#########################################################################

ulimit -s unlimited #Without this line of code you will get segmentation fault trying to run lga

#Remember that the pdbs need to be in the same folder as the collect script (collect_PDB.pl)

# $1 is predicted_pdb_file_location from the python file GDT_TS_Automation line 72, (e.g. ranked_0, ranked_1 ... unrelaxed_5)
# $2 is template_pdb_file_location from the python file GDT_TS_Automation line 72, (e.g. 6t1z)

#First, we make a copy of the files (the result from AlphaFold and the Template) into the MOL2 folder
#$(basename $1) is to get only the last part of the path from $1 which is predicted_pdb_file_location
#$(basename $2) is to get only the last part of the path from $1 which is template_pdb_file_location
cp $1 LGA_Zemla/MOL2/$(basename $1)
cp $2 LGA_Zemla/MOL2/$(basename $2)


#Due to how Zemla's code works we need to go to the directory, instead of calling it from the base script
#Otherwise, it will send an error that the file is not found
cd LGA_Zemla/MOL2
./collect_PDB.pl $(basename $1) > pdb1_plus_pdb2.pdb
./collect_PDB.pl $(basename $2) >> pdb1_plus_pdb2.pdb

#Now we go to the folder in which the 'lga' binary is
cd ..
#Remember that ./lga reads the file directly from the MOL2 folder, so we do not need to specify the full path, just the filename
#So in this case an example of call would be:
#-4 -o2 -gdc -lga_m -stral -d:4.0 pdb1_plus_pdb2
#'pdb1_plus_pdb2' being the result of calling 'collect_PDB_pl' two times to get a single file with both molecules
#This single file is what we feed to the lga algorithm.
####################### IMPORTANT ########################
#####										                             #####
#####  pdb1_plus_pdb2 needs to be in the MOL2 folder #####
#####										                             #####
##########################################################
./lga -4 -o2 -gdc -lga_m -stral -d:4.0 pdb1_plus_pdb2.pdb > MOL2/superposition_$(basename $1)



#In order to run the GDT_TS calculation which is done with the following parameters:
#./lga -3 -o2 -gdc -lga_m -stral -d:4.0 -al [superposition_file_name]
#The superposition file alone is not enough, we need to curate it by aggregating the information from pdb1_plus_pdb2
#So, we read the information from 'pdb1_plus_pdb2' into a variable and then feed it to the file
#Thus, we rewrite the file with the code below
#Where the variable pdb1_plus_pdb2 is assigned the value of the file pdb1_plus_pdb2 ( yes both have the same names dont get confused )
pdb1_plus_pdb2=$(<MOL2/pdb1_plus_pdb2.pdb) 
#We open the file 'MOL2/superposition_$(basename $1)' and append the information of the variable pdb1_plus_pdb2
echo "$pdb1_plus_pdb2" >> MOL2/superposition_$(basename $1)

#Now the file is ready for the last part of the process, the GDT_TS score
#We call the lga binary file below with the newly curated 'superposition_$(basename $1)' file
#Remember that ./lga reads the file directly from the MOL2 folder, so we do not need to specify the full path, just the filename
#In this case it would be 'superposition_$(basename $1)'
./lga -3 -o2 -gdc -lga_m -stral -d:4.0 -al superposition_$(basename $1) >> MOL2/GDT_TS_$(basename $1)


