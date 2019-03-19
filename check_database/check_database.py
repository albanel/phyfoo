from ModuleFonction_BuildOrthocis import *
import os

original_generic_tables = ["current_genome_db","current_parameter","current_member","current_method_link","current_method_link_species_set","current_species_set","current_homology"]

def check_generic_tables(host,database,user,pw,db,bHaveMySQLdb,checkdir) :
	generic_dico = {}
	for table in original_generic_tables :
		query = "DESCRIBE "+table+" ;"
        	table_description = SubmitMySQLCommand(query,host,database,user,pw,db,bHaveMySQLdb,False)
		fields = []
        	for tRow in table_description :
                	field = tRow[0]
			fields.append(field)
			generic_dico[table] = fields
	print "++++++++++ DATABASE GENERIC FIELDS  ++++++++++"
	print generic_dico
	download_file = checkdir+"/check_database.sh"
	os.system("sh "+ download_file + " " + checkdir)
	with open(checkdir+"/generic_dico.txt","r") as generic_dico_file : 
		print "++++++++++ ENSEMBL GENERIC FIELDS  ++++++++++"
		for line in generic_dico_file :
			print line
