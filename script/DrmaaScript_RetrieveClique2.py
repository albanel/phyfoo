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
iIndexLine=int(sys.argv[2]) #line number for agruments to use
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
		sListOfComponentId=tLine[3:] #table for all first Gene Id of component connex research

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

db=None
if bHaveMySQLdb:
	db = MySQLdb.connect(host=host,user=user,passwd=pw,db=database) 
########################################################################

selectfieldid="gene_id"
if "global" in sGlobalTag:
                selectfieldid="global_id"

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
	retrieve_clique_bash_file_name = "quick_clique_for_cc_"+str(sComponentId)+".sh"
	retrieve_clique_bash_file = open(retrieve_clique_bash_file_name,"w")
	cmd = quick_clique_dir+"/bin/qc --algorithm=degeneracy --input-file=$1 > ${1}.cl"
	retrieve_clique_bash_file.write(cmd)
	retrieve_clique_bash_file.close()
	print "Calculating cliques..."
	os.system("bash "+retrieve_clique_bash_file_name+" "+output_file_name)
	clique_file_name = output_file_name+".cl"
	clique_file = open(clique_file_name,"r")
	decoded_clique_file_name = "decoded_"+clique_file_name
	decoded_clique_file = open(decoded_clique_file_name,"w")
	clique_file.readline()
	clique_file.readline()
	clique_ct = 0
	print "**********quick-cliques**********"
	oCliqueList = []
	for clique in clique_file :
		clique_ct += 1
		clique = clique.split(" ")
		clique_list = []
		for num in clique :
			num = int(num)
			gene_id = dico_id_to_node[num]
			clique_list.append(gene_id)
		#print sorted(clique_list)
		oCliqueList.append(sorted(clique_list))
		#decoded_clique_file.write(str(clique_list)+"\n")
	clique_file.close()	
	decoded_clique_file.close()
	print "nb cliques : "+str(clique_ct)
	'''
	print "**********networkx**********"
	#Construct Graph
	Graph=nx.Graph()
	Graph.add_nodes_from(tComponentContent)
	Graph.add_edges_from(tLinkOfComponent)
	#Research clique and add to tables
	oCliqueList = nx.find_cliques(Graph)
	nb_clique = 0
	oCliqueList = sorted(oCliqueList)
	'''
	for oClique in oCliqueList:
		print sorted(oClique)
		iCompteurClique+=1	
		#print "new clique "+str(iCompteurClique)
		#sCliqueId=sComponentId+"_"+str(iCompteurClique)
		tPositivGenomeDbId=[]
		tEquivalentTableBool=[]
		tEquivalentTableElmt=[]
		tTaxonList=[]
		for sEntityId in oClique:
			print("A4%s" % (sEntityId))
			try:
				#print "\t"+sEntityId
				sGenomeDbId=dLinkEntityId_specieId[str(sEntityId)]
				print("A42")
				sCodeId="sp"+sGenomeDbId
				tPositivGenomeDbId.append(sCodeId)
				print("A5")
				tEquivalentTableBool.append("1")
				tEquivalentTableElmt.append(sEntityId)
				tTaxonList.append(sGenomeDbId)
			except :
			 	print  sys.exc_info()	
				sys.exit("WARNING:gene_id possibly not found in table gene for gene_id (or global_id) ="+str(sEntityId))
			
		sTargetedColumn=",".join(tPositivGenomeDbId)
		#sBooleanValue=",".join(tEquivalentTableBool)
		sElmtValue=",".join(tEquivalentTableElmt)
	 
		print("A6")
		#sSQLaddNewLine="insert into "+sBdDVersion+"_clique_bool ("+sTargetedColumn+") values ("+sBooleanValue+");"
		#SubmitMySQLCommand(sSQLaddNewLine,host,database,user,pw,db,bHaveMySQLdb,False)

		sSQLaddNewLine2="insert into "+sBdDVersion+sGlobalTag+"_clique_element ("+sTargetedColumn+") values ("+sElmtValue+");"
		trySubmitMySQLCommand(sSQLaddNewLine2,host,database,user,pw,db,bHaveMySQLdb,False)

FILE=open("End.txt","w")
FILE.close()

