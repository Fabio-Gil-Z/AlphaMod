import os
import glob
import time
import argparse
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import NoSuchElementException
from selenium.common.exceptions import StaleElementReferenceException
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.support import expected_conditions
from Screenshot import Screenshot_Clipping

# QMEAN: https://swissmodel.expasy.org/qmean/
# ProsaWeb: https://prosa.services.came.sbg.ac.at/prosa.php
# MolProbity: http://molprobity.biochem.duke.edu/
# Procheck: https://saves.mbi.ucla.edu/

parser = argparse.ArgumentParser(description="Using web automation with Selenium")

parser.add_argument('-Q', '--QMEAN', metavar='', required=False, help="If True, executes the QMEAN automation")
parser.add_argument('-P', '--PROSA', metavar='', required=False, help="If True, executes the PROSA automation")
parser.add_argument('-M', '--MOLPROBITY', metavar='', required=False, help="If True, executes the QMEAN automation")
parser.add_argument('-C', '--PROCHECK', metavar='', required=False, help="If True, executes the PROCHECK automation")
parser.add_argument('-T', '--TARGETNAME', metavar='', required=False, help="Target Name")

args = parser.parse_args()

def Qmean(pdb_file_path):
	cwd = os.getcwd()
	Google_Chrome_Selenium_Web_Bot_Folder = f'{cwd}/Structure_Assessment_Section/Google_Chrome_Selenium_Web_Bot'
	basename = os.path.basename(pdb_file_path)
	message = f"Calculating QMEAN {os.path.basename(args.TARGETNAME)}"
	lenght_basename = len(basename)
	lenght_message = len(message)
	spaces = 3
	full_message = "#"*spaces + " " + message + " " + basename + " " + "#"*spaces
	print(f"#"*spaces + "#"*lenght_message + "#"*(lenght_basename+3) +"#"*spaces)
	print(f"#"*spaces + " "*lenght_message + " "*(lenght_basename+3) +"#"*spaces)
	print(full_message)
	print(f"#"*spaces + " "*lenght_message + " "*(lenght_basename+3) +"#"*spaces)
	print(f"#"*spaces + "#"*lenght_message + "#"*(lenght_basename+3) +"#"*spaces)

	ignored_exceptions=(NoSuchElementException, StaleElementReferenceException, TimeoutException)
	driver = webdriver.Chrome(f'{Google_Chrome_Selenium_Web_Bot_Folder}/chromedriver')  # Optional argument, if not specified will search path.
	driver.get('https://swissmodel.expasy.org/qmean/')

	# Configuring the driver to start maximized
	driver.maximize_window()

	# The 'select_coordinate_file_button' can be found at the website https://swissmodel.expasy.org/qmean/
	driver.find_element(By.NAME, "structureFile").clear()
	select_coordinate_file_button = driver.find_element(By.NAME, "structureFile")


	# Send the file location to the button (the PDB fle)
	select_coordinate_file_button.send_keys(pdb_file_path)

	# Submit button
	submit_button = driver.find_element(By.ID, "submitButton")
	# submit_button = driver.find_element_by_id('submitButton')
	time.sleep(10) # wait for the file to load
	submit_button.click()

	# Wait until the job is done
	# WebDriverWait(driver, 600).until(EC.text_to_be_present_in_element((By.TAG_NAME,"body"),"QMEANDisCo Global:"))
	WebDriverWait(driver, 600, ignored_exceptions=ignored_exceptions)\
                        .until(expected_conditions.text_to_be_present_in_element((By.TAG_NAME,"body"),"QMEANDisCo Global:"))
	superposition_window = driver.window_handles[-1]
	driver.switch_to.window(superposition_window)

	QMEAN_score_raw_website_text = driver.find_element(By.TAG_NAME, "body").text.split(sep='\n')

	for line in QMEAN_score_raw_website_text:
		if "QMEANDisCo Global:" in line:
			QMEAN_score = line
			if not os.path.exists(f'{Google_Chrome_Selenium_Web_Bot_Folder}/QMEAN'):
				os.mkdir(f'{Google_Chrome_Selenium_Web_Bot_Folder}/QMEAN')
				print("\nDirectory: ", f'{Google_Chrome_Selenium_Web_Bot_Folder}/QMEAN', " has been Created ")

			QMEAN_score = QMEAN_score.split(sep=":")
			QMEAN_score = QMEAN_score[1]
			QMEAN_score = QMEAN_score.split(sep=" ")
			QMEAN_score = QMEAN_score[0]
			# cleaned .csv file
			if not os.path.isfile(f'{Google_Chrome_Selenium_Web_Bot_Folder}/QMEAN/QMEAN_results_cleaned_csv.csv'):
				# print("FILE DOES NOT EXIST I WRITE")
				with open(f'{Google_Chrome_Selenium_Web_Bot_Folder}/QMEAN/QMEAN_results_cleaned_csv.csv', 'w') as file:
					file.write(f'Model,Q-MEAN_DisCo_Global\n')
					file.write(f'{os.path.basename(pdb_file_path)},{QMEAN_score}\n')
			else:
				# print("FILE EXIST I APPEND")
				with open(f'{Google_Chrome_Selenium_Web_Bot_Folder}/QMEAN/QMEAN_results_cleaned_csv.csv', 'a') as file:
					file.write(f'{os.path.basename(pdb_file_path)},{QMEAN_score}\n')

			# Define the path to save the screenshot
			path_to_save_body_screenshot = f'{Google_Chrome_Selenium_Web_Bot_Folder}/QMEAN'

			#Saving screenshot
			save_screenshot(driver, pdb_file_path, QMEAN_score, path_to_save_body_screenshot)

	# time.sleep(5) # Let the user actually see something!
	driver.quit() # quit and close

def PROSA_Z_SCORE(pdb_file_path):
	cwd = os.getcwd()
	Google_Chrome_Selenium_Web_Bot_Folder = f'{cwd}/Structure_Assessment_Section/Google_Chrome_Selenium_Web_Bot'
	basename = os.path.basename(pdb_file_path)
	message = f"Calculating PROSA {os.path.basename(args.TARGETNAME)}"
	lenght_basename = len(basename)
	lenght_message = len(message)
	spaces = 3
	full_message = "#"*spaces + " " + message + " " + basename + " " + "#"*spaces
	print(f"#"*spaces + "#"*lenght_message + "#"*(lenght_basename+3) +"#"*spaces)
	print(f"#"*spaces + " "*lenght_message + " "*(lenght_basename+3) +"#"*spaces)
	print(full_message)
	print(f"#"*spaces + " "*lenght_message + " "*(lenght_basename+3) +"#"*spaces)
	print(f"#"*spaces + "#"*lenght_message + "#"*(lenght_basename+3) +"#"*spaces)

	driver = webdriver.Chrome(f'{Google_Chrome_Selenium_Web_Bot_Folder}/chromedriver')  # Optional argument, if not specified will search path.
	driver.get('https://prosa.services.came.sbg.ac.at/prosa.php')

	# Configuring the driver to start maximized
	driver.maximize_window()

	# The 'userfile' can be found at the website https://prosa.services.came.sbg.ac.at/prosa.php
	driver.find_element(By.NAME, "userfile").clear()
	select_userfile_button = driver.find_element(By.NAME, "userfile")


	# Send the file location to the button (the PDB fle)
	select_userfile_button.send_keys(pdb_file_path)

	# Submit button
	submit_button = driver.find_element(By.NAME, "submButton")
	# submit_button = driver.find_element_by_id('submitButton')
	time.sleep(5) # wait for the file to load
	submit_button.click()

	# Wait until the job is done
	WebDriverWait(driver, 600).until(EC.text_to_be_present_in_element((By.TAG_NAME,"body"),"Z-Score:"))
	superposition_window = driver.window_handles[-1]
	driver.switch_to.window(superposition_window)

	prosa_z_score_raw_website_text = driver.find_element(By.TAG_NAME, "body").text.split(sep='\n')

	# print(prosa_z_score_raw_website_text)
	# for i in prosa_z_score_raw_website_text:
	# 	print(i)
	# print(type(prosa_z_score_raw_website_text))


	for line in prosa_z_score_raw_website_text:
		if "Z-Score:" in line:
			PROSA = line
			if not os.path.exists(f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROSA'):
				os.mkdir(f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROSA')
				print("\nDirectory: ", f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROSA', " has been Created ")

			PROSA = PROSA.split(sep=":")
			PROSA = PROSA[1].replace(" ", "")
			# raw .csv file
			if not os.path.isfile(f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROSA/PROSA_results_raw_csv.csv'):
				# print("FILE DOES NOT EXIST I WRITE")
				with open(f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROSA/PROSA_results_raw_csv.csv', 'w') as file:
					file.write(f'Model,PROSA\n')
					file.write(f'{os.path.basename(pdb_file_path)},{PROSA}\n')
			else:
				# print("FILE EXIST I APPEND")
				with open(f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROSA/PROSA_results_raw_csv.csv', 'a') as file:
					file.write(f'{os.path.basename(pdb_file_path)},{PROSA}\n')


			# Define the path to save the screenshot
			path_to_save_body_screenshot = f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROSA'

			#Saving screenshot
			save_screenshot(driver, pdb_file_path, PROSA, path_to_save_body_screenshot)



	# time.sleep(5) # Let the user actually see something!
	driver.quit() # quit and close


def MOLPROBITY(pdb_file_path):
	# chrome_options = webdriver.ChromeOptions()
	# chrome_options.add_experimental_option("detach", True)
	# driver = webdriver.Chrome(f'{os.getcwd()}/chromedriver', options=chrome_options)  # Optional argument, if not specified will search path.
	# driver.get('http://molprobity.biochem.duke.edu/')

	cwd = os.getcwd()
	Google_Chrome_Selenium_Web_Bot_Folder = f'{cwd}/Structure_Assessment_Section/Google_Chrome_Selenium_Web_Bot'
	basename = os.path.basename(pdb_file_path)
	message = f"Calculating MOLPROBITY {os.path.basename(args.TARGETNAME)}"
	lenght_basename = len(basename)
	lenght_message = len(message)
	spaces = 3
	full_message = "#"*spaces + " " + message + " " + basename + " " + "#"*spaces
	print(f"#"*spaces + "#"*lenght_message + "#"*(lenght_basename+3) +"#"*spaces)
	print(f"#"*spaces + " "*lenght_message + " "*(lenght_basename+3) +"#"*spaces)
	print(full_message)
	print(f"#"*spaces + " "*lenght_message + " "*(lenght_basename+3) +"#"*spaces)
	print(f"#"*spaces + "#"*lenght_message + "#"*(lenght_basename+3) +"#"*spaces)

	driver = webdriver.Chrome(f'{Google_Chrome_Selenium_Web_Bot_Folder}/chromedriver')  # Optional argument, if not specified will search path.
	driver.get('http://molprobity.biochem.duke.edu/')

	# Configuring the driver to start maximized
	driver.maximize_window()

	# The 'uploadFile' can be found at the website http://molprobity.biochem.duke.edu/
	driver.find_element(By.NAME, "uploadFile").clear()
	select_uploadFile_button = driver.find_element(By.NAME, "uploadFile")


	# Send the file location to the button (the PDB fle)
	select_uploadFile_button.send_keys(pdb_file_path)

	# Submit button
	submit_button = driver.find_element(By.NAME, "cmd")
	# submit_button = driver.find_element_by_id('submitButton')
	time.sleep(3) # wait for the file to load
	submit_button.click()

	# Wait until the job is done
	WebDriverWait(driver, 600).until(EC.text_to_be_present_in_element((By.TAG_NAME,"body"),"Uploaded PDB file as"))
	superposition_window = driver.window_handles[-1]
	driver.switch_to.window(superposition_window)

	# Getting the 'Continue >' button
	select_continue_button = driver.find_element(By.NAME, "cmd")
	select_continue_button.click()

	# Wait until the job is done
	WebDriverWait(driver, 600).until(EC.text_to_be_present_in_element((By.TAG_NAME,"body"),"Analyze geometry without all-atom contacts"))
	superposition_window = driver.window_handles[-1]
	driver.switch_to.window(superposition_window)


	# Getting the 'Analyze geometry' link
	driver.find_element("link text", "Analyze geometry without all-atom contacts").click()

	# Wait until the job is done
	WebDriverWait(driver, 600).until(EC.text_to_be_present_in_element((By.TAG_NAME,"body"),"Analyze all-atom contacts and geometry"))
	superposition_window = driver.window_handles[-1]
	driver.switch_to.window(superposition_window)

	################################## ---- CLASHES & CLASSCORE ---- ##################################
	# Getting the 'Run programs to perform these analyses >' button
	select_clash_score_checkbox = driver.find_element(By.NAME, "chartClashlist")
	select_clash_score_checkbox.click()

	################################## ---- CLASHES & CLASSCORE ---- ##################################

	# Getting the 'Run programs to perform these analyses >' button
	select_run_programs_button = driver.find_element(By.NAME, "cmd")
	select_run_programs_button.click()

	alert = driver.switch_to.alert
	alert.accept()

	WebDriverWait(driver, 600).until(EC.text_to_be_present_in_element((By.TAG_NAME,"body"),"Single-criterion visualizations"))
	superposition_window = driver.window_handles[-1]
	driver.switch_to.window(superposition_window)

	# Values found at the end of the query, we are getting the colors
	# = driver.find_element_by_xpath("").value_of_css_property('background-color')

	color_list = []

	MOLPROBITY_website_text = driver.page_source.split(sep='\n')



	#color_count = getColor(MOLPROBITY_website_text)
	#print(color_count)

	molprobity_score = getMolprobityScore(MOLPROBITY_website_text)
	molprobity_score = molprobity_score.split()
	molprobity_score = molprobity_score[2]

	if not os.path.exists(f'{Google_Chrome_Selenium_Web_Bot_Folder}/MOLPROBITY'):
		os.mkdir(f'{Google_Chrome_Selenium_Web_Bot_Folder}/MOLPROBITY')
		print("\nDirectory: ", f'{Google_Chrome_Selenium_Web_Bot_Folder}/MOLPROBITY', " has been Created ")

	# raw .csv file
	if not os.path.isfile(f'{Google_Chrome_Selenium_Web_Bot_Folder}/MOLPROBITY/MOLPROBITY_results_raw_csv.csv'):
		# print("FILE DOES NOT EXIST I WRITE")
		with open(f'{Google_Chrome_Selenium_Web_Bot_Folder}/MOLPROBITY/MOLPROBITY_results_raw_csv.csv', 'w') as file:
			file.write(f'Model,MOLPROBITY\n')
			file.write(f'{os.path.basename(pdb_file_path)},{molprobity_score}\n')
	else:
		# print("FILE EXIST I APPEND")
		with open(f'{Google_Chrome_Selenium_Web_Bot_Folder}/MOLPROBITY/MOLPROBITY_results_raw_csv.csv', 'a') as file:
			file.write(f'{os.path.basename(pdb_file_path)},{molprobity_score}\n')


	# Define the path to save the screenshot
	path_to_save_body_screenshot = f'{Google_Chrome_Selenium_Web_Bot_Folder}/MOLPROBITY'

	#Saving screenshot
	save_screenshot(driver, pdb_file_path, {molprobity_score}, path_to_save_body_screenshot)

	driver.quit() # quit and close


def PROCHECK(pdb_file_path):
	cwd = os.getcwd()
	Google_Chrome_Selenium_Web_Bot_Folder = f'{cwd}/Structure_Assessment_Section/Google_Chrome_Selenium_Web_Bot'
	basename = os.path.basename(pdb_file_path)
	message = f"Calculating PROCHECK {os.path.basename(args.TARGETNAME)}"
	lenght_basename = len(basename)
	lenght_message = len(message)
	spaces = 3
	full_message = "#"*spaces + " " + message + " " + basename + " " + "#"*spaces
	print(f"#"*spaces + "#"*lenght_message + "#"*(lenght_basename+3) +"#"*spaces)
	print(f"#"*spaces + " "*lenght_message + " "*(lenght_basename+3) +"#"*spaces)
	print(full_message)
	print(f"#"*spaces + " "*lenght_message + " "*(lenght_basename+3) +"#"*spaces)
	print(f"#"*spaces + "#"*lenght_message + "#"*(lenght_basename+3) +"#"*spaces)

	driver = webdriver.Chrome(f'{Google_Chrome_Selenium_Web_Bot_Folder}/chromedriver')  # Optional argument, if not specified will search path.
	driver.get('https://saves.mbi.ucla.edu/')

	# Configuring the driver to start maximized
	driver.maximize_window()

	# The 'select_coordinate_file_button' can be found at the website https://swissmodel.expasy.org/qmean/
	driver.find_element(By.NAME, "pdbfile").clear()
	select_choose_file_button = driver.find_element(By.NAME, "pdbfile")

	# Send the file location to the button (the PDB fle)
	select_choose_file_button.send_keys(pdb_file_path)

	# Submit button
	submit_button = driver.find_element(By.ID, "startjob")
	# submit_button = driver.find_element_by_id('submitButton')
	time.sleep(3) # wait for the file to load
	submit_button.click()

	# Wait until the job is done

	WebDriverWait(driver, 600).until(EC.presence_of_element_located((By.ID, "procheck")))
	# WebDriverWait(driver, 600).until(EC.text_to_be_present_in_element((By.TAG_NAME,"body"),"References"))
	superposition_window = driver.window_handles[-1]
	driver.switch_to.window(superposition_window)


	# Start button
	start_button = driver.find_element(By.ID, "procheck")
	# start_button = driver.find_element_by_id('startButton')
	time.sleep(3) # wait for the file to load
	start_button.click()

	# Wait until the job is done
	WebDriverWait(driver, 600).until(EC.text_to_be_present_in_element((By.TAG_NAME,"body"),"Results"))
	superposition_window = driver.window_handles[-1]
	driver.switch_to.window(superposition_window)

	# //*[@id="wprocheck"]/div/center/span/a
	ignored_exceptions=(NoSuchElementException, StaleElementReferenceException, TimeoutException)
	# results_button = driver.find_element_by_xpath(r'//*[@id="wprocheck"]/div/center/span/a')


	WebDriverWait(driver, 600, ignored_exceptions=ignored_exceptions)\
                        .until(expected_conditions.presence_of_element_located((By.XPATH, "//*[text()='Results']")))


	while True:
		try:
			results_button = driver.find_element(By.XPATH, "//*[text()='Results']").click()
			break
		except:
			pass


	superposition_window = driver.window_handles[-1]
	driver.switch_to.window(superposition_window)
	procheck_raw_website_text = driver.find_element(By.TAG_NAME, "body").text.split(sep='\n')
	#print(procheck_raw_website_text)

	# for i in procheck_raw_website_text:
	# 	print(i)

	for line in procheck_raw_website_text:
		if "Ramachandran plot:" in line:
			PROCHECK = line
			PROCHECK = PROCHECK.replace('+', "")
			PROCHECK = PROCHECK.replace('|', "")
			PROCHECK = PROCHECK.replace('*', "")
			PROCHECK = PROCHECK.split()
			PROCHECK = f'{PROCHECK[2][:-1]}|{PROCHECK[4][:-1]}|{PROCHECK[6][:-1]}|{PROCHECK[8][:-1]}'

			if not os.path.exists(f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROCHECK'):
				os.mkdir(f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROCHECK')
				print("Directory: ", f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROCHECK', " has been Created ")

			# raw .csv file
			if not os.path.isfile(f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROCHECK/PROCHECK_results_raw_csv.csv'):
				# print("FILE DOES NOT EXIST I WRITE")
				with open(f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROCHECK/PROCHECK_results_raw_csv.csv', 'w') as file:
					file.write(f'Model,PROCHECK\n')
					file.write(f'{os.path.basename(pdb_file_path)},{PROCHECK}\n')
			else:
				# print("FILE EXIST I APPEND")
				with open(f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROCHECK/PROCHECK_results_raw_csv.csv', 'a') as file:
					file.write(f'{os.path.basename(pdb_file_path)},{PROCHECK}\n')

			# Define the path to save the screenshot
			path_to_save_body_screenshot = f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROCHECK'

			# Saving screenshot, had to use built in screenshot because scrolling screnshot was not working
			# The images were being chop off at the right side
			output_file_name = os.path.basename(pdb_file_path)
			output_file_name = os.path.splitext(output_file_name)[0]
			driver.save_screenshot(f'{path_to_save_body_screenshot}/{output_file_name} {PROCHECK}')

			#save_screenshot(driver, pdb_file_path, {PROCHECK}, path_to_save_body_screenshot)
	driver.quit() # quit and close


def getMolprobityScore(html_source_code_as_list):
	for html_string in html_source_code_as_list:
				if "MolProbity score" in html_string:
					index_before_score = html_string.find("\">")
					score = html_string[index_before_score+2:index_before_score+6]
					return f'MolProbity score {score}'


def getColor(html_source_code_as_list):
	HEXADECIMAL_COLOR_LENGTH = 7

	RED_COLOR = "#ff9999"
	YELLOW_COLOR = "#ffff99"
	GREEN_COLOR = "#99ff99"

	RED = 0
	YELLOW = 0
	GREEN = 0

	list_of_parameters = ["Poor rotamers",
	"Favored rotamers",
	"Ramachandran outliers",
	"Ramachandran favored",
	"Rama distribution Z-score",
	"deviations",
	"Bad bonds",
	"Bad angles",
	"Cis Prolines",
	"Twisted Peptides",
	"CaBLAM outliers",
	"CA Geometry outliers",
	"Chiral volume outliers"]

	for parameter in list_of_parameters:
		for html_string in html_source_code_as_list:
			if parameter in html_string:
				color = html_string[html_string.find("#"):html_string.find("#")+HEXADECIMAL_COLOR_LENGTH]
				if color == RED_COLOR:
					RED += 1
				elif color == YELLOW_COLOR:
					YELLOW += 1
				elif color == GREEN_COLOR:
					GREEN += 1
	return f'RED: {RED}, YELLOW: {YELLOW}, GREEN: {GREEN}'



def save_screenshot(driver, pdb_file_path, score, path_to_save_body_screenshot):
	#Saving screenshot
	output_file_name = os.path.basename(pdb_file_path)
	output_file_name = os.path.splitext(output_file_name)[0]
	ob=Screenshot_Clipping.Screenshot()
	img=ob.full_Screenshot(driver,save_path=path_to_save_body_screenshot,image_name=f'{output_file_name}_{score}.png')

def main():
	cwd = os.getcwd()
	Google_Chrome_Selenium_Web_Bot_Folder = f'{cwd}/Structure_Assessment_Section/Google_Chrome_Selenium_Web_Bot'
	pdb_file_path_list = glob.glob(f'{Google_Chrome_Selenium_Web_Bot_Folder}/pdb_files/*.pdb')
	pdb_file_path_list.sort()

	#########--- tmp_ouputfiles ---#########
	'''
    tmp_ouputfiles ---> if len(tmp_ouputfiles) == 6
    Means the Web Crawler was able to get all the information
    5 Screenshots from 5 pdbs and 1 csv file containing the scores
    If for some reason the Web Crawler crashes, (Sever side, TimeOutException, Stale or any other inconvenience)
    The while true will keep iterating until it is able to get all the 6 required files
    '''
	
	# QMEAN
	if args.QMEAN.lower() == "true":
		while True:
			for pdb_file_path in pdb_file_path_list:
				print("\n\nFile location:")
				print(f'{pdb_file_path}\n')
				Qmean(pdb_file_path)
			tmp_ouputfiles = glob.glob(f'{Google_Chrome_Selenium_Web_Bot_Folder}/QMEAN/*')
			if len(tmp_ouputfiles) == 6:
				break
	# PROSA
	if args.PROSA.lower() == "true":
		while True:
			for pdb_file_path in pdb_file_path_list:
				print("\n\nFile location:")
				print(f'{pdb_file_path}\n')
				PROSA_Z_SCORE(pdb_file_path)
			tmp_ouputfiles = glob.glob(f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROSA/*')
			if len(tmp_ouputfiles) == 6:
				break
	
	# MOLPROBITY
	if args.MOLPROBITY.lower() == "true":
		while True:
			for pdb_file_path in pdb_file_path_list:
				print("\n\nFile location:")
				print(f'{pdb_file_path}\n')
				MOLPROBITY(pdb_file_path)
			tmp_ouputfiles = glob.glob(f'{Google_Chrome_Selenium_Web_Bot_Folder}/MOLPROBITY/*')
			if len(tmp_ouputfiles) == 6:
				break
	
	# PROCHECK
	if args.PROCHECK.lower() == "true":
		while True:
			for pdb_file_path in pdb_file_path_list:
				print("\n\nFile location:")
				print(f'{pdb_file_path}\n')
				PROCHECK(pdb_file_path)
			tmp_ouputfiles = glob.glob(f'{Google_Chrome_Selenium_Web_Bot_Folder}/PROCHECK/*')
			if len(tmp_ouputfiles) == 6:
				break

if __name__=='__main__':
	main()
