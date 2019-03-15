import os
from ModuleFonction_BuildOrthocis import *

def check_database(host,database,user,pw,db,bHaveMySQLdb,False,genericdir) :
	#step 1 : count the columns for all tables 
	dTableToColumnNumber={}
        query = "SHOW TABLES"
        tables = SubmitMySQLCommand(query,host,database,user,pw,db,bHaveMySQLdb,True)
        for tRow in tables :
                table = tRow[0]
		query = "SELECT COUNT(*) FROM INFORMATION_SCHEMA.COLUMNS WHERE table_name = '"+table+"' ;"
		nber_columns = SubmitMySQLCommand(query,host,database,user,pw,db,bHaveMySQLdb,True)
		for nbre in nber_columns :
			dTableToColumnNumber[table]=nbre[0]
	print dTableToColumnNumber
	#step 2 : count fields of corresponding files
	#step 2.1 : count fields of generic tables (Ensembl Compara)
	generic_files_list = os.listdir(genericdir)
	dFileToFields={}
	for generic_file in generic_files_list :
		if ".txt.gz" in generic_file :
			unziped_name = generic_file.replace(".gz","")
			if unziped_name not in generic_files_list :
				dzip_cmd = "zcat "+genericdir+"/"+generic_file+" > "+genericdir+"/"+unziped_name
				os.system(dzip_cmd)
			with open(genericdir+"/"+unziped_name,"r") as current_generic_file :
				for line in current_generic_file :
					line = line.split("\t")
					print unziped_name
					nbre_fields=len(line)
					print nbre_fields
					break
			generic_file = generic_file.replace(".txt.gz","")
			generic_file = generic_file.replace("95","current")
			dFileToFields[generic_file]=nbre_fields
	print dFileToFields
	#step3 : compare both dictionnary
	for key in dFileToFields :
		if key in dTableToColumnNumber.keys() :
			if dFileToFields[key] == dTableToColumnNumber[key] :
				print "OK for file and table "+key
			else :
				print "Warning ! Difference between file "+key+" and its table " 
	
		else :
			print "The file " + key + " doesn't correspond to any table in the database"
