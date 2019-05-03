import sys
import os
import networkx as nx

try:
	import MySQLdb
	bHaveMySQLdb=True
except ImportError:
	bHaveMySQLdb=False
	
from ModuleFonction_BuildOrthocis import *
#force mysql command
bHaveMySQLdb=False
#Arguments
sPathToArgFile=sys.argv[1] #path to the file of arguments
iIndexLine=int(sys.argv[2]) # number for agruments to use
########################################################################
iCurrentIndex=0
for sLine in open(sPathToArgFile):
	if iIndexLine==iCurrentIndex:
		sLine=sLine.replace("\n","")
		tLine=sLine.split("\t")

		if len(tLine)==0:
			exit("No Component Id for search component connex")

		sBdDVersion=tLine[0] #internal prefixe of the bdd
		sPathToConfFile=tLine[1] #path to the conf file
		sGlobalTag=tLine[2] #add globalTag to table name
		SInsert=tLine[3] #compute or insertion, if false only compute, if true only insertion
		sListOfComponentId=tLine[4:] #table for all first Gene Id of component connex research

		break
	iCurrentIndex+=1
########################################################################
#Mysql
conffile=sPathToConfFile
oConfig =populateGeneralConf(conffile)
host  = oConfig.get('cf_host')
database = oConfig.get('cf_database')
user = oConfig.get('cf_user')
pw = oConfig.get('cf_pw')

quick_clique_dir = oConfig.get('cf_quick_clique_dir')

clique_job_id = iIndexLine+1

db=None
if bHaveMySQLdb:
	db = MySQLdb.connect(host=host,user=user,passwd=pw,db=database) 
########################################################################

selectfieldid="gene_id"
if "global" in sGlobalTag:
                selectfieldid="global_id"

if "ok_for_compute" in os.listdir(".") :
	print "Cliques already computed for the job_id "+str(clique_job_id)

if "ok_for_insertion" in os.listdir(".") :
	print "Cliques already written in the database for the job_id "+str(clique_job_id)

#SInsert = "False"

if SInsert == "False" and "ok_for_compute" not in os.listdir('.') :
	print "Insert = False. Only compute cliques."
	decoded_clique_file_name = "total_cliques_for_cc"+str(clique_job_id)
	decoded_clique_file = open(decoded_clique_file_name,"w")
	ct = 0
	for sComponentId in sListOfComponentId :	
		ct += 1
		print "Retriving links for connex component "+str(sComponentId)
		file_name = "links_for_cc_"+str(sComponentId)+".txt"
		linkFile = open(file_name,"w")
		tComponentContent=[]
		tLinkOfComponent=[]
		dLinkEntityId_specieId={}
		iCompteurClique=0
		allMb=""
		#Get all geneId of the component
		sSQLRetrieveAllGene="select "+selectfieldid+" from "+sBdDVersion+sGlobalTag+"_component_connex where component_id="+sComponentId+";"
		tResult=trySubmitMySQLCommand(sSQLRetrieveAllGene,host,database,user,pw,db,bHaveMySQLdb,True)
		for tRow in tResult:
			tComponentContent.append(tRow[0])
		if len(tComponentContent)>0:	
			allMb="("+",".join(tComponentContent)+")"
		else :
			sys.exit("void cc for id "+sComponentId)
		#print(tComponentContent,allMb)	
		#Retrieve genomeDbId of each gene
		sSQLRetrieveGenomeDbId="select "+selectfieldid+",genome_db_id from "+sBdDVersion+"_gene where "+selectfieldid+" in "+allMb+";"
		#print sSQLRetrieveGenomeDbId
		tResult=trySubmitMySQLCommand(sSQLRetrieveGenomeDbId,host,database,user,pw,db,bHaveMySQLdb,True)
		#print tResult
		for Row in tResult:
			dLinkEntityId_specieId[Row[0]]=Row[1]
			#print Row[0]+"->"+dLinkEntityId_specieId[Row[0]]
		#print dLinkEntityId_specieId
		#Get all link between gene
		sSQLRetrieveAllInterestHomolink="select * from "+sBdDVersion+sGlobalTag+"_homolink where  first in "+allMb+" and second in "+allMb+"  ;"
		#print  sSQLRetrieveAllInterestHomolink	
		tResult=trySubmitMySQLCommand(sSQLRetrieveAllInterestHomolink,host,database,user,pw,db,bHaveMySQLdb,True)
		for tRow in tResult:
			tLinkOfComponent.append((tRow[0],tRow[1]))
			linkFile.write(tRow[0]+"\t"+tRow[1]+"\n")
		linkFile.close()
		
		liste_elt = []

		input_file_name = file_name

		with open(input_file_name,"r") as input_file :
			nbre_edges = 0
			for line in input_file :
				line = line.replace("\n","")
				line = line.split("\t")
				liste_elt.append(line[0])
				liste_elt.append(line[1])
				nbre_edges += 1

		liste_elt = set(liste_elt)

		print("Nbre de vertex : "+str(len(liste_elt)))

		dico_node_to_id = {}
		dico_id_to_node = {}

		for i, item in enumerate(liste_elt) :
			dico_node_to_id[item] = i
			dico_id_to_node[i] = item


		output_file_name = input_file_name.replace(".txt","")+"_treated.txt"

		print("Creating "+output_file_name+"...")

		output_file = open(output_file_name,"w")

		with open(input_file_name,"r") as input_file :
			output_file.write(str(len(liste_elt))+"\n")
			output_file.write(str(nbre_edges*2)+"\n")
			for line in input_file :
				line = line.replace("\n","")
				line = line.split("\t")
				output_file.write(str(dico_node_to_id[line[0]])+","+str(dico_node_to_id[line[1]])+"\n")
				output_file.write(str(dico_node_to_id[line[1]])+","+str(dico_node_to_id[line[0]])+"\n")

		output_file.close()
		print "Calculating cliques..."
		print "**********quick-cliques**********"
		retrieve_clique_bash_file_name = "quick_clique_for_cc_"+str(sComponentId)+".sh"
		retrieve_clique_bash_file = open(retrieve_clique_bash_file_name,"w")
		cmd = quick_clique_dir+"/bin/qc --algorithm=degeneracy --input-file=$1 > ${1}.cl"
		retrieve_clique_bash_file.write(cmd)
		retrieve_clique_bash_file.close()
		os.system("bash "+retrieve_clique_bash_file_name+" "+output_file_name)
		clique_file_name = output_file_name+".cl"
		clique_file = open(clique_file_name,"r")
		clique_file.readline()
		clique_file.readline()
		clique_ct = 0
		oCliqueList = []
		print "Decoding clique file..."
		for clique in clique_file :
			clique_ct += 1
			clique = clique.split(" ")
			clique_list = []
			for num in clique :
				num = int(num)
				gene_id = dico_id_to_node[num]
				clique_list.append(gene_id)
			#print sorted(clique_list)
			#oCliqueList.append(sorted(clique_list))
			current_line = ""
			for gene in sorted(clique_list) :
				current_line += str(gene)+","
			current_line = current_line.rstrip(",")
			decoded_clique_file.write(current_line+"\n")
		clique_file.close()	
		#print "nb cliques : "+str(clique_ct)
		#Clean directory
		print "Cleaning up directory..."
		os.remove(file_name)
		os.remove(output_file_name)
		os.remove(clique_file_name)
		os.remove(retrieve_clique_bash_file_name)
		'''
		print "**********networkx**********"
		#Construct Graph
		Graph=nx.Graph()
		Graph.add_nodes_from(tComponentContent)
		Graph.add_edges_from(tLinkOfComponent)
		#Research clique and add to tables
		oCliqueList2 = nx.find_cliques(Graph)
		for clique in oCliqueList2 :
			print(sorted(clique))
		'''
	decoded_clique_file.close()
	ok_for_compute = open("ok_for_compute","w")
	ok_for_compute.close()

#SInsert = "True"

if SInsert == "True" and "ok_for_insertion" not in os.listdir('.') :
	sSQL_delete_clique_for_current_job = "DELETE FROM "+sBdDVersion+sGlobalTag+"_clique_element WHERE clique_job_id = "+str(clique_job_id)+" ;"
	SubmitMySQLCommand(sSQL_delete_clique_for_current_job,host,database,user,pw,db,bHaveMySQLdb,True)
	print "Insert = True. Clique files will be treated and inserted in the database."
	if "End.txt" in os.listdir('.') :
		os.remove("End.txt")
	target_file = "total_cliques_for_cc"+str(clique_job_id)
	if target_file not in os.listdir('.') :
		sys.exit("no clique file to treat for id"+str(clique_job_id))
	else :
		print "Creating table file for cc"+str(clique_job_id)
		sSQLretrieve_species_columns = "SELECT column_name FROM INFORMATION_SCHEMA.COLUMNS WHERE table_name = 'current_global_clique_element';"
		species_columns = SubmitMySQLCommand(sSQLretrieve_species_columns,host,database,user,pw,db,bHaveMySQLdb,True)

		clique_file = open(target_file,"r")
		reconsitued_table_file_name = "table_for_job_"+str(clique_job_id)
		reconsitued_table_file = open(reconsitued_table_file_name,"w")
		clique_ct = 0
		for oClique in clique_file :
			clique_ct += 1
			oClique = oClique.replace("\n","")
			oClique = oClique.split(",")
			specie_to_gene = {}
			line_to_insert = str(clique_job_id)
			for gene in oClique :
				specie_id = gene[0:3]
				specie_to_gene[specie_id] = gene
			for column in species_columns :
				if column[0] != "clique_job_id" :
					specie_id = column[0].replace('sp','')
					try :
						line_to_insert += ","+specie_to_gene[specie_id]
					except KeyError :
						line_to_insert += ",NULL"
			#print(line_to_insert)
			reconsitued_table_file.write(line_to_insert+"\n")
			if clique_ct%100000 == 0 :
				print "Treated cliques : "+str(clique_ct)
		reconsitued_table_file.close()
		clique_file.close()
	#sSQL_load_data = "LOAD DATA LOCAL INFILE '"+reconsitued_table_file_name+"' INTO TABLE " +sBdDVersion+sGlobalTag+"_clique_element FIELDS TERMINATED BY ',' LINES TERMINATED BY '\\n';"
	#SubmitMySQLCommand(sSQL_load_data,host,database,user,pw,db,bHaveMySQLdb,True)
	ok_for_insertion = open("ok_for_insertion","w")
	ok_for_insertion.close()

FILE=open("End.txt","w")
FILE.close()
