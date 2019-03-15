import sys
import shutil
import os
import re
import cPickle

from optparse import OptionParser


from ModuleFTP import *
from ModuleCreateTable_BuildOrthocis import *
from ModuleFonction_BuildOrthocis import *

#######################################



#Options
parser = OptionParser()
parser.add_option("-s","--step", dest="step")
parser.add_option("-c","--conf", dest="conf")
parser.add_option("-l","--logdir", dest="ldirectory")
(options, args) = parser.parse_args()

#O2M:TODO:0.1
bUSE_ONE_2_MANY=True

logDir="log"
if options.ldirectory is not None:
	logDir=options.ldirectory


base=os.path.basename(options.step)
sStdOutFile=logDir+"/"+os.path.splitext(base)[0]+"_info.log"


conffile=options.conf
if not conffile:
        print("Error : no configuration file -c  defined, process broken")
        sys.exit()

################################
#recover file
sDebugFile=options.step
bDebugFile=False
sCurrentDebugFile="RecoveryFile.step"
if sDebugFile:
        bDebugFile=True

        tDebugFile=sDebugFile.split("/")
        if len(tDebugFile)>1:
                #sPath="/".join(tDebugFile[:-1])
		sPath=logDir
                sFileName=tDebugFile[-1]
                sCurrentDebugFile=sPath+"/outRecov"+sFileName
        else:
                sCurrentDebugFile="outRecov"+sDebugFile

DebugFile=open(sCurrentDebugFile,"w")
DebugFile.write("#DebugFile to set the variables controlling tasks")
DebugFile.close()
################################

try:
        os.makedirs(logDir)
	print ("dir created,  logs available in "+logDir)
except OSError :
        print (" logs available in "+logDir)

bPrint=True
dotruDRMAA=False


if os.path.exists(sStdOutFile):
	os.remove(sStdOutFile)
FILE=open(sStdOutFile,"w")
FILE.close()

try:
	import MySQLdb
	bHaveMySQLdb=True
	printDebug("----------Check if module MySQLdb is present : Yes",sStdOutFile,bPrint)
except ImportError:
	bHaveMySQLdb=False
	printDebug("----------Check if module MySQLdb is present : No",sStdOutFile,bPrint)

#force mysl command usage : 
bHaveMySQLdb=False


bHaveDrmaa=False
if dotruDRMAA:
	#Already false ???
	try:
		import drmaa
		bHaveDrmaa=True
		printDebug("----------Check if module Drmma is present : Yes",sStdOutFile,bPrint)
	except ImportError:
		printDebug("----------Check if module Drmma is present : No",sStdOutFile,bPrint)
	

########################################################################
#SETUP
sDefaultValueForMissignDebugOption="true"

bDoRequestEnsembl=True #Ok
if bDoRequestEnsembl:
	bDoDownload=True #Ok
	bDoRequestForData=True #Ok
	bDoRequestForSpecies=True #Ok
else:
	bDoDownload=bDoRequestForData=bDoRequestForSpecies=False

bDoRequestHGNC=True #Ok


bDoAllbase=True
if bDoAllbase:
	bDoGenericTable=True #Ok
	bGunzipGenomeDb=True #Ok
	bFillTableGenomeDbTable=True #Ok
	bDoSpeciesSpecificTable=True #Ok 
	bFillAutoGenerateTable=True #Ok
	if bFillAutoGenerateTable:
		bFillToolTable=True
		bFillGenericTable=True #Ok
		bFillSpecificTable=True #Ok 
	else:
		bFillToolTable=bFillGenericTable=bFillDeductibleTable=False
	bFillHomolinkTable=True 
else:
	bDoGenericTable=bGunzipGenomeDb=bFillTableGenomeDbTable=bDoSpeciesSpecificTable=bFillAutoGenerateTable=bFillDeductibleTable=False
	
bRetrieveComponentConnex=True
if bRetrieveComponentConnex:
	bCreateComponentConnexTable=True #Ok
	bFillComponentConnexTable=True #Todo : Test
else:
	bCreateComponentConnexTable=bFillComponentConnexTable=False

bRetrieveClique=True
if bRetrieveClique:
	bCreateCliqueBoolTable=True #Ok
	bFillCliqueTables=True #Todo : Test
else:
	bCreateCliqueBoolTable=bFillCliqueTables=False

bBuildHGNCtables=True
if bBuildHGNCtables:
	bBuildHGNC=True #Ok
	bBuildHGNCassoClique=True #Ok
else:
	bBuildHGNC=bBuildHGNCassoClique=False
	
bCreateFasta=True #Todo : Test
########################################################################
#Options

#sDebugFile=options.step
#bDebugFile=False
#sCurrentDebugFile="RecoveryFile.step"
#
#if sDebugFile:
#        bDebugFile=True
#
  #      tDebugFile=sDebugFile.split("/")
 #       if len(tDebugFile)>1:
   #             sPath="/".join(tDebugFile[:-1])
#                sFileName=tDebugFile[-1]
#                sCurrentDebugFile=sPath+"/outRevoc"+sFileName
#        else:
 #               sCurrentDebugFile="outRevoc"+sDebugFile

#DebugFile=open(sCurrentDebugFile,"w")
#DebugFile.write("#DebugFile to set the variables controlling tasks")
#DebugFile.close()

conffile=options.conf
conffileabs=os.path.abspath(conffile)

if not conffile:
        printDebug("Error : no configuration file -c  defined, process broken",sStdOutFile,bPrint)
        sys.exit()

# read variables from config file
###mysql
oConfig =populateGeneralConf(conffile)


host  = oConfig.get('cf_host')
database = oConfig.get('cf_database')
user = oConfig.get('cf_user')
pw = oConfig.get('cf_pw')

speciesdata = oConfig.get('cf_species2data')

species2data=readJsonstring(speciesdata)



upstreamMode=oConfig.get('cf_upstreamMode')
if not upstreamMode:
	upstreamMode=1
printDebug("upstreamMode::"+str(upstreamMode),sStdOutFile,bPrint)
biotype=oConfig.get('cf_biotype')
if not biotype:
	biotype='protein_coding'
printDebug("biotype::"+str(biotype),sStdOutFile,bPrint)

###remote files
ftpurl=oConfig.get('cf_ftpurl')
ftpHGNC=oConfig.get('cf_ftpHGNC')
HGNCcompleteSetPath=oConfig.get('cf_HGNCcompleteSetPath')
sFolderForPersistentData=oConfig.get('cf_sFolderForPersistentData')
sQueue=oConfig.get('cf_sQueue')
sQueueCC=oConfig.get('cf_sQueueCC')
sQueueClique=oConfig.get('cf_sQueueClique')
jobexecdir=oConfig.get('cf_jobexecdir')
###########application related parameters
iMaxJobForParallelize=int(oConfig.get('cf_maxJob'))-1
iMaxJobForParallelizeHLink=int(oConfig.get('cf_maxJobHLink'))-1
excludedGID=oConfig.get('cf_excludedGID')
#excludedGIDAr=[]
#if excludedGID and len(excludedGID)>1:
#		for e in excludedGID.split(","):
#			excludedGIDAr.append(e.strip())
		
#Execute the process on all or particular species
speciesfilter=oConfig.get('cf_speciesfilter')
if not speciesfilter:
        speciesfilter=".*"
#Ensembl version
sVersionEnsembl=oConfig.get('cf_version')
sVersionOrthocis="current"
sVersionTemp=sVersionOrthocis

if not sVersionEnsembl:
        printDebug("Error : no version in conf file  , process broken",sStdOutFile,bPrint)
        sys.exit()

#threshold (size of upstream dna sequence to extracted
threshold=oConfig.get('cf_threshold')
threshold=int(threshold)
###locale work and download directories

oudir=oConfig.get('cf_workdir')
scriptdir=oConfig.get('cf_scriptdir')

if not os.path.exists(oudir):
	try:    
		os.makedirs(oudir)
	except OSError :
		printDebug("can not create "+oudir+" ",sStdOutFile,bPrint)

scompco=scriptdir+"/"+"compco.pl"
shomolinkA=scriptdir+"/"+"homolinkA.py"
shomolinkB=scriptdir+"/"+"homolinkB.py"
shomolinkC=scriptdir+"/"+"homolinkC.py"
sdefineGlobalID=scriptdir+"/"+"defineGlobalID.py"


sDrmaaScript_ExtractFromFasta=scriptdir+"/"+"ExtractFromFastaPARA.py"
sDrmaaScript_RetrieveComponentConnex=scriptdir+"/"+"DrmaaScript_RetrieveComponentConnex.py"
sDrmaaScript_RetrieveClique=scriptdir+"/"+"DrmaaScript_RetrieveClique.py"
#dRelationGenomeDBId_SpecieName="Empty"

sPathToHGNCfolder=oudir+"/HGNCdata"
ensdir=oudir+"/ens/"
persisDir=oudir+"/"+sFolderForPersistentData+"/persist/"
fastaInDir=oudir+"/"+sFolderForPersistentData+"/fastain/"
fastaOutDir=oudir+"/"+sFolderForPersistentData+"/fastaout/"

persistdictoFile = persisDir+"/RelationGenomeDbId_SpecieName.dic"

for edir in [oudir+"/"+sFolderForPersistentData, persisDir, fastaInDir,fastaOutDir,ensdir]:
	if not os.path.exists(edir):
	 	try:
			os.makedirs(edir)
		except OSError :
			printDebug("can not create "+edir+" ",sStdOutFile,bPrint)


########################################################################
#Use personnal DebugFile
tOptionList=["DoRequestEnsembl","DoDownload","DoRequestForData","DoRequestForSpecies",
"DoRequestHGNC",
"DoAllbase","DoGenericTable","GunzipGenomeDb","FillTableGenomeDbTable","DoSpeciesSpecificTable","FillAutoGenerateTable","FillGenericTable","FillSpecificTable","FillHomolinkTable",
"RetrieveComponentConnex","CreateComponentConnexTable","FillComponentConnexTable",
"RetrieveClique","CreateCliqueBoolTable","FillCliqueTables",
"BuildHGNCtables","BuildHGNC","BuildHGNCassoClique",
"CreateFasta"
]
if bDebugFile:
	printDebug("----------Load DebugFile option",sStdOutFile,bPrint)
	oConfig = CustomConfig(sDebugFile)
	
	for sOption in tOptionList:
		try:
			sValue=oConfig.get(sOption)
		except ConfigParser.NoOptionError:
			sValue=sDefaultValueForMissignDebugOption
		bBooleanValue=TranslateInBool(sValue)
		if bBooleanValue==TypeError:
			exit("DebugFile : option "+sOption+", unknow value : "+sValue+"\nNeed a boolean")
			sys.exit()
		sNameofBooleanVariable="b"+sOption
		locals()[sNameofBooleanVariable]=bBooleanValue

		#print sNameofBooleanVariable,":",locals()[sNameofBooleanVariable]

########################################################################

#mysql connection

db=None

########################################################################
#Conf for Test
tCommandSQL_Test=[
"CREATE TABLE "+sVersionTemp+"""_clique_member (
  `clique_id` int(10) unsigned NOT NULL,
  `gene_id` int(10) unsigned NOT NULL,
  KEY `clique_id_idx` (`clique_id`),
  KEY `gene_id_idx` (`gene_id`)
)ENGINE=MyISAM DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+sVersionTemp+"""_clique ;
CREATE TABLE  $nameCliquetemp  (clique_id  int(10) unsigned NOT NULL AUTO_INCREMENT,
   taxon_list  varchar(1000), PRIMARY KEY (clique_id) ) 
ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+sVersionTemp+"""_gene;
CREATE TABLE `$nameGenetemp` (
  `gene_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `stable_id` varchar(128) NOT NULL,
  `genome_db_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`gene_id`),
  KEY (`stable_id`) 
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1 MAX_ROWS=100000000;"""
]
tListeGeneTest=[
519897,
297694,
528194,
509690,
223193,
232030,
377171,
443229,
452519,
501508,
19600,
351685
]
########################################################################
#Main
'''
sgenomeIDNSQLCmd="SELECT genome_db_id,name FROM "+sVersionOrthocis+"_genome_db "
sgenomeIDNSQLCmd+=" WHERE taxon_id IS NOT NULL "
if excludedGID and len(excludedGID)>1:
	sgenomeIDNSQLCmd+="and genome_db_id not in ("+excludedGID+")"
sgenomeIDNSQLCmd+=";"

tSpeciesTable=""
try:
  # fm 01 2017
  # select species set  for all steps
  tSpeciesTable=SubmitMySQLCommand(sgenomeIDNSQLCmd,host,database,user,pw,db,bHaveMySQLdb,True)


  dRelationGenomeDBId_SpecieName={}
  for tRow in tSpeciesTable:
	iGenome_db_id=tRow[0]
	sName=tRow[1]
	#if iGenome_db_id not in excludedGIDAr:
	dRelationGenomeDBId_SpecieName[sName]=int(iGenome_db_id)
        printDebug("====" +sName+":"+iGenome_db_id+" ",sStdOutFile,bPrint)
  printDebug("-Persist genomeDbId/SpecieName relation : Starting",sStdOutFile,bPrint)
  cPickle.dump(dRelationGenomeDBId_SpecieName,open(persistdictoFile,"w"))
  printDebug("-Persist genomeDbId/SpecieName : Ending",sStdOutFile,bPrint)

except:
        printDebug("table not created yet ! (normal for step1 ) ",sStdOutFile,bPrint)
'''

def dico():
	sgenomeIDNSQLCmd="SELECT genome_db_id,name FROM "+sVersionOrthocis+"_genome_db "
	sgenomeIDNSQLCmd+=" WHERE taxon_id IS NOT NULL "
	if excludedGID and len(excludedGID)>1:
        	sgenomeIDNSQLCmd+="and genome_db_id not in ("+excludedGID+")"
	sgenomeIDNSQLCmd+=";"
	tSpeciesTable=""
	try:
  	# fm 01 2017
  	# select species set  for all steps
  		tSpeciesTable=SubmitMySQLCommand(sgenomeIDNSQLCmd,host,database,user,pw,db,bHaveMySQLdb,True)
  		dRelationGenomeDBId_SpecieName={}
  		for tRow in tSpeciesTable:
        		iGenome_db_id=tRow[0]
        		sName=tRow[1]
        		#if iGenome_db_id not in excludedGIDAr:
        		dRelationGenomeDBId_SpecieName[sName]=int(iGenome_db_id)
        		printDebug("====" +sName+":"+iGenome_db_id+" ",sStdOutFile,bPrint)
  		printDebug("-Persist genomeDbId/SpecieName relation : Starting",sStdOutFile,bPrint)
  		cPickle.dump(dRelationGenomeDBId_SpecieName,open(persistdictoFile,"w"))
  		printDebug("-Persist genomeDbId/SpecieName : Ending",sStdOutFile,bPrint)
		global tSpeciesTable
        	global dRelationGenomeDBId_SpecieName
		global sgenomeIDNSQLCmd
	except:
        	printDebug("table not created yet ! (normal for step1 ) ",sStdOutFile,bPrint)
dico()	
##albane Jai deplace la partie precedente car le dictionnaire doit etre cree plus tard sinon souci lors de la creation des tables specifiques
 #Retrieve all files from Ensembl
if bDoRequestEnsembl:
	printDebug("----------Step 1 RequestEnsembl : Starting",sStdOutFile,bPrint)
	
	if bDoRequestForData:
		printDebug("- Download Global file : Starting",sStdOutFile,bPrint)
		log=ftpSpeciesRelationDownload(bDoDownload,ftpurl,ensdir,speciesfilter,sVersionEnsembl)
		printDebug("["+",".join(str(oLog) for oLog in log)+"]",sStdOutFile,bPrint)
		printDebug("- Download Global file : Ending",sStdOutFile,bPrint)
		#DEBUG---------------------------------------------------
		DEBUG="""
		ls [ensdir]
			ensembl_compara_[sVersionEnsembl]
		ls [ensdir]/ensembl_compara_[sVersionEnsembl]
			[sVersionEnsembl]_genome_db.txt.gz	   [sVersionEnsembl]_method_link_species_set.txt.gz
			[sVersionEnsembl]_homology_member.txt.gz   [sVersionEnsembl]_method_link.txt.gz
			[sVersionEnsembl]_homology.txt.gz	   [sVersionEnsembl]_species_set.txt.gz
			[sVersionEnsembl]_gene_member.txt.gz
		"""
		#DEBUG---------------------------------------------------
	else:
		printDebug("- Download Global file : No",sStdOutFile,bPrint)

	printDebug("DoRequestForData=false",sCurrentDebugFile,False)

	if bDoRequestForSpecies:
		printDebug("- FASTA here --- Download Species file : Starting",sStdOutFile,bPrint)
		log=ftpSpeciesDataDownload(bDoDownload,ftpurl,ensdir,fastaInDir,speciesfilter,sVersionEnsembl)
		printDebug("["+",".join(str(oLog) for oLog in log)+"]",sStdOutFile,bPrint)
		printDebug("- Download Species file : Ending",sStdOutFile,bPrint)
		#DEBUG---------------------------------------------------
		DEBUG="""
		ls [ensdir]
			ensembl_compara_[sVersionEnsembl]
			[specie1]_core_[sVersionEnsembl]_[int]	   [specie2]_core_[sVersionEnsembl]_[int]
			...
		ls [ensdir]/[specie1]_core_[sVersionEnsembl]_[int]
			seq_region.txt.gz	    transcript.txt.gz
			assembly.txt.gz  	    coord_system.txt.gz
			[specie1].[shortname].[sVersionEnsembl].dna.toplevel.fa.gz

			WARNING : toplevel and/or primary_assembly fasta.gz file possible. At least one must be present
		"""
		#DEBUG---------------------------------------------------
	else:	

		printDebug("- Download Species file : No",sStdOutFile,bPrint)

	printDebug("DoRequestForSpecies=false",sCurrentDebugFile,False)

	printDebug("----------Step 1 RequestEnsembl : Ending",sStdOutFile,bPrint)
else:
	printDebug("----------Step 1 RequestEnsembl : No",sStdOutFile,bPrint)

printDebug("DoRequestEnsembl=false",sCurrentDebugFile,False)

if bDoRequestHGNC:
	printDebug("----------Step 2 RequestHGNC : Starting",sStdOutFile,bPrint)
	log=ftpHGNCfileDownload(bDoDownload,ftpHGNC,sPathToHGNCfolder,HGNCcompleteSetPath)
	printDebug("["+",".join(str(oLog) for oLog in log)+"]",sStdOutFile,bPrint)
	printDebug("----------Step 2 RequestHGNC : Ending",sStdOutFile,bPrint)
else:
	printDebug("----------Step 2 RequestHGNC : No",sStdOutFile,bPrint)
	printDebug("DoRequestHGNC=false",sCurrentDebugFile,False)

#Create Allbase
if bDoAllbase:

	#Create generic table ie not need the genome_db of specie in her name
	if bDoGenericTable:
		printDebug("----------Step 3 Build generic table : Starting",sStdOutFile,bPrint)
		#O2M:TODO:0.4
		tCommandSQL=DefineSQLRequestGeneric(sVersionOrthocis,str(threshold))
		for sSQLcmd in tCommandSQL:
			SubmitMySQLCommand(sSQLcmd,host,database,user,pw,db,bHaveMySQLdb,False)
		printDebug("----------Step 3 Build generic table : Ending",sStdOutFile,bPrint)
		#DEBUG---------------------------------------------------
		DEBUG="""
		mysql --batch -u [user] -p[pw] -h [host] [database] -e "show tables;"
			[sVersionOrthocis]_HGNC
			[sVersionOrthocis]_asso_HGNC_clique
			[sVersionOrthocis]_clique
			[sVersionOrthocis]_clique_member
			[sVersionOrthocis]_exec_param
			[sVersionOrthocis]_gene
			[sVersionOrthocis]_genome_db
			[sVersionOrthocis]_homolink
			[sVersionOrthocis]_homology
			[sVersionOrthocis]_homology_member
			[sVersionOrthocis]_job
			[sVersionOrthocis]_member
			[sVersionOrthocis]_method_link
			[sVersionOrthocis]_method_link_species_set
			[sVersionOrthocis]_motif
			[sVersionOrthocis]_parameter
			[sVersionOrthocis]_scan
			[sVersionOrthocis]_sequence
			[sVersionOrthocis]_species_set
			[sVersionOrthocis]_tool
			[sVersionOrthocis]_type
		"""
		#DEBUG---------------------------------------------------
	else:
		printDebug("----------Step 3 Build generic table : No",sStdOutFile,bPrint)

	printDebug("DoGenericTable=false",sCurrentDebugFile,False)
	
	if bGunzipGenomeDb:
		printDebug("----------Step 4 Unzip GenomeDb File : starting",sStdOutFile,bPrint)
		#Unzip Ensembl file and rename genome_db for easy target with regex
		GetCommandGunzip_genomeDB(ensdir,sVersionEnsembl)

		printDebug("----------Step 4 Unzip GenomeDb File : Ending",sStdOutFile,bPrint)
		#DEBUG---------------------------------------------------
		DEBUG="""
		ls [ensdir]/ensembl_compara_[sVersionEnsembl]
			[sVersionEnsembl]_genome_db.bdd	   [sVersionEnsembl]_method_link_species_set.txt.gz
			[sVersionEnsembl]_homology_member.txt.gz  [sVersionEnsembl]_method_link.txt.gz
			[sVersionEnsembl]_homology.txt.gz	   [sVersionEnsembl]_species_set.txt.gz
			[sVersionEnsembl]_gene_member.txt.gz
		"""
		#DEBUG---------------------------------------------------
	else:
		printDebug("----------Step 4 Unzip GenomeDb File : No",sStdOutFile,bPrint)

	printDebug("GunzipGenomeDb=false",sCurrentDebugFile,False)

	#Fill genome_db table for obtain the relation species name/genome_db id
	if bFillTableGenomeDbTable:
		printDebug("----------Step 5 Fill GenomeDb table : Starting",sStdOutFile,bPrint)
		sEnsemblComparaFolder=ensdir+"/ensembl_compara_"+sVersionEnsembl
		ExecuteMySQLImportCommand(host,database,user,pw,sEnsemblComparaFolder+"/"+"genome_db.txt",sVersionOrthocis+"_genome_db")

		printDebug("----------Step 5 Fill GenomeDb table : Ending",sStdOutFile,bPrint)
		#DEBUG---------------------------------------------------
		DEBUG="""
		mysql -u [user] -p[pw] -h [host] [database] -e "select * from [sVersionOrthocis]_genome_db limit 2;"
			+--------------+----------+-------------------+----------+------------------+-----------------+---------+
			| genome_db_id | taxon_id | name              | assembly | assembly_default | genebuild       | locator |
			+--------------+----------+-------------------+----------+------------------+-----------------+---------+
			|            4 |    31033 | takifugu_rubripes | FUGU4    |                1 | 2007-11-Ensembl | 0       |
			|           27 |    51511 | ciona_savignyi    | CSAV2.0  |                1 | 2006-04-Ensembl | 0       |
			+--------------+----------+-------------------+----------+------------------+-----------------+---------+
		"""
		#DEBUG---------------------------------------------------
	else:

		printDebug("----------Step 5 Fill GenomeDb table : No",sStdOutFile,bPrint)

	printDebug("FillTableGenomeDbTable=false",sCurrentDebugFile,False)
	#Create species table ie table who need the genome_db id in her name
	if bDoSpeciesSpecificTable:
		dico()
		#fm 01 2017
		#tSQLtable=SubmitMySQLCommand(sgenomeIDNSQLCmd,host,database,user,pw,db,bHaveMySQLdb,True)
		liste_espece=speciesfilter.split(",") 
		for tRow in tSpeciesTable:
			iGenome_db_id=tRow[0]
			sName=tRow[1]
			#if iGenome_db_id not in excludedGIDAr:
			#dRelationGenomeDBId_SpecieName[sName]=int(iGenome_db_id)
			'''
			p = re.compile(speciesfilter)
			if p.match(sName):
				 tSpecificRequestForSpecies=DefineSQLRequestSpecificToSpecie(sVersionOrthocis,str(iGenome_db_id),str(threshold) )
				 #idx=0;
				 for sSQLcommandPart2 in tSpecificRequestForSpecies:
					##########					
					#idx+=1
					#if(idx==1):
					#	print sSQLcommandPart2
					#	SubmitMySQLCommand(sSQLcommandPart2,host,database,user,pw,db,bHaveMySQLdb,False)
					SubmitMySQLCommand(sSQLcommandPart2,host,database,user,pw,db,bHaveMySQLdb,False)
			'''                          
			for espece in liste_espece :
			        p = re.compile(espece)
                        	if p.match(sName):
                                	 tSpecificRequestForSpecies=DefineSQLRequestSpecificToSpecie(sVersionOrthocis,str(iGenome_db_id),str(threshold))
                                	 #idx=0;
                                	 for sSQLcommandPart2 in tSpecificRequestForSpecies:
                                       	 ##########                                      
                                       	 #idx+=1
                                       	 #if(idx==1):
                                       	 #       print sSQLcommandPart2
                                       	 #       SubmitMySQLCommand(sSQLcommandPart2,host,database,user,pw,db,bHaveMySQLdb,False)
                                         	SubmitMySQLCommand(sSQLcommandPart2,host,database,user,pw,db,bHaveMySQLdb,False)

		printDebug("----------Step 6 Build species-named table : Ending",sStdOutFile,bPrint)
		#printDebug("-Persist genomeDbId/SpecieName relation : Starting",sStdOutFile,bPrint)
		#cPickle.dump(dRelationGenomeDBId_SpecieName,open(persistdictoFile,"w"))
		#printDebug("-Persist genomeDbId/SpecieName : Ending",sStdOutFile,bPrint)
		#DEBUG---------------------------------------------------
		DEBUG="""
		mysql -u [user] -p[pw] -h [host] [database] -e "show tables;"
			+---------------------------------+
			| Tables_in_vapordev              |
			+---------------------------------+
			| ...       	                  |
			| current_assembly_109            |
			| ...   	                  |
			| current_coord_system_109        |
			| current_dna_109                 |
			| ...              	          |
			| current_gene_109                |
			| ...                  	          |
			| current_seq_region_109          |
			| ...         	                  |
			| current_transcript_109          |
			| ...	                          |
			+---------------------------------+
		"""
		#DEBUG---------------------------------------------------
	else:

		printDebug("----------Step 6 Build species-named table : No",sStdOutFile,bPrint)
		#printDebug("-Load genomeDbId/SpecieName relation : Starting",sStdOutFile,bPrint)
		#try:
		#	dRelationGenomeDBId_SpecieName=cPickle.load(open(persistdictoFile,"r"))
		#	printDebug("-Load genomeDbId/SpecieName relation : Ending",sStdOutFile,bPrint)
		#except:
		#	exit("\nImpossible to load : "+persistdictoFile+"\nProcess broken")
		#	sys.exit()

	printDebug("DoSpeciesSpecificTable=false",sCurrentDebugFile,False)

	if bFillAutoGenerateTable:
		printDebug("----------Step 7 Fill auto-generate table : Starting",sStdOutFile,bPrint)
		if bFillToolTable:
			printDebug("-Fill tool table : Starting",sStdOutFile,bPrint)
			dTool2Id=GetContentOfToolTable()
			for sTool in dTool2Id.keys():
				sSQLFillToolTable="insert ignore into "+sVersionOrthocis+"_tool values("+str(dTool2Id[sTool][0])+",'"+sTool+"','"+dTool2Id[sTool][1]+"');"
				SubmitMySQLCommand(sSQLFillToolTable,host,database,user,pw,db,bHaveMySQLdb,False)
			printDebug("-Fill tool table : Ending",sStdOutFile,bPrint)
		if bFillGenericTable:
			#first level : generic file (homology_member, etc)
			printDebug("-Fill generic table : Starting",sStdOutFile,bPrint)
			#O2M:TODO:0.2
			ExecuteMultiMySQLImportFromAFolder(host,database,user,pw,ensdir+"/ensembl_compara_"+sVersionEnsembl,sVersionOrthocis,sVersionEnsembl,bUSE_ONE_2_MANY)

			printDebug("-Fill generic table : Ending",sStdOutFile,bPrint)
			#DEBUG---------------------------------------------------
			DEBUG="""
			mysql -u [user] -p[pw] -h [host] [database] -e "select * from [sVersionOrthocis]_homology_member limit 2;"
				+-------------+-----------+-------------------+------------+----------+---------+----------+
				| homology_id | member_id | peptide_member_id | cigar_line | perc_cov | perc_id | perc_pos |
				+-------------+-----------+-------------------+------------+----------+---------+----------+
				|   300000001 | 300000001 |         300000002 | 207M       |      100 |     100 |      100 |
				|   300000001 | 300000003 |         300000004 | 207M       |      100 |     100 |      100 |
				+-------------+-----------+-------------------+------------+----------+---------+----------+
			"""
			#DEBUG---------------------------------------------------

		else:
			printDebug("-Fill generic table : No",sStdOutFile,bPrint)

		printDebug("FillGenericTable=false",sCurrentDebugFile,False)
		 
		if bFillSpecificTable:
			#second level : for each specie, go in folder and treat specific species file
			printDebug("-Fill specific table : Starting",sStdOutFile,bPrint)
			tFileList=os.listdir(ensdir)
			for specieFolder in tFileList:
				if "core" in specieFolder:
					sSpeciePath=ensdir+"/"+specieFolder
					sNameSpecie=specieFolder.split("_core")[0]
					if sNameSpecie in dRelationGenomeDBId_SpecieName.keys():					
						iGenomeDbId=dRelationGenomeDBId_SpecieName[sNameSpecie]
					        printDebug("- file "+sSpeciePath+" as dict entry "+sNameSpecie+" "+str(iGenomeDbId),sStdOutFile,bPrint)	
						ExecuteMultiMySQLImportFromAFolder(host,database,user,pw,sSpeciePath,sVersionOrthocis,sVersionEnsembl,bUSE_ONE_2_MANY,str(iGenomeDbId))				
					else:
						printDebug("-warning file "+sSpeciePath+" as no corresponding dict entry "+sNameSpecie,sStdOutFile,bPrint)
				
			printDebug("-Fill specific table : Ending",sStdOutFile,bPrint)
		else:
			printDebug("-Fill specific table : No",sStdOutFile,bPrint)

		printDebug("FillSpecificTable=false",sCurrentDebugFile,False)

		printDebug("----------Step 7 Fill auto-generate table : Ending",sStdOutFile,bPrint)
		#DEBUG---------------------------------------------------
		DEBUG="""
		mysql -u [user] -p[pw] -h [host] [database] -e "select * from [sVersionOrthocis]_seq_region_109 limit 2;"
			+---------------+------------+-----------------+--------+
			| seq_region_id | name       | coord_system_id | length |
			+---------------+------------+-----------------+--------+
			|             1 | GL193649.1 |               1 | 568412 |
			|             2 | GL193838.1 |               1 | 461840 |
			+---------------+------------+-----------------+--------+
		"""
		#DEBUG---------------------------------------------------
	else:
		printDebug("----------Step 7 Fill auto-generate table : No",sStdOutFile,bPrint)

	printDebug("FillAutoGenerateTable=false",sCurrentDebugFile,False)

	if bFillHomolinkTable:
		#TAG1	
		#fm 03 2016
		
		dico()    
		printDebug("----------Step 8.00 filter homology_member : Starting",sStdOutFile,bPrint)
		sqmhm="DROP TABLE IF EXISTS "+sVersionOrthocis+"_tmp_hm ;"
                SubmitMySQLCommand(sqmhm,host,database,user,pw,db,bHaveMySQLdb,False)
                sqmhm="create table "+sVersionOrthocis+"_tmp_hm  select * from "+sVersionOrthocis+"_homology_member;"
                SubmitMySQLCommand(sqmhm,host,database,user,pw,db,bHaveMySQLdb,False) 
		sqmhm="delete from "+sVersionOrthocis+"_homology_member;"
		SubmitMySQLCommand(sqmhm,host,database,user,pw,db,bHaveMySQLdb,False)
                sqmhm="insert into "+sVersionOrthocis+"_homology_member select  m.*  from "+sVersionOrthocis+"_tmp_hm m, "+sVersionOrthocis+"_homology h  where m.homology_id =h.homology_id ;"
               	SubmitMySQLCommand(sqmhm,host,database,user,pw,db,bHaveMySQLdb,False) 
		sqmhm="drop table "+sVersionOrthocis+"_tmp_hm;"
                SubmitMySQLCommand(sqmhm,host,database,user,pw,db,bHaveMySQLdb,False)
		printDebug("----------Step 8.00 filter homology_member : end",sStdOutFile,bPrint)		
	
		printDebug("----------Step 8 Fill Homolink table : Starting",sStdOutFile,bPrint)
		SQLgetOrthologMethodLinkId="select method_link_id  from "+sVersionOrthocis+"_method_link where type='ENSEMBL_ORTHOLOGUES';"
		sMethodId=None
		tSQLtable=SubmitMySQLCommand(SQLgetOrthologMethodLinkId,host,database,user,pw,db,bHaveMySQLdb,True)
		for uniqRow in tSQLtable:
			sMethodId=uniqRow[0]
			printDebug("----------Step 8 --0--",sStdOutFile,bPrint)

		tListOfGenomeDbId=sorted(dRelationGenomeDBId_SpecieName.values()) 

		tArrayOfArg=[]
		tListOfCouple=[]
 		printDebug("----------Step 8 --1-- : Ending",sStdOutFile,bPrint)
###############"

#Add AAL
		sTempName1=getTempdir()+"TempArgFile"+str(random.random())+str(random.random())+".txt"
		sTempName2=getTempdir()+"TempArgFile"+str(random.random())+str(random.random())+".txt"
		iLineNumber1=0
		iLineNumber2=0
		printDebug("----------Step 8 --2-- : Ending",sStdOutFile,bPrint)
		ArgFile1=open(sTempName1,"w")

		for iIndexGenome1 in range(0,len(tListOfGenomeDbId)-1):
			sGenomeDbId1=tListOfGenomeDbId[iIndexGenome1]
			for iIndexGenome2 in range(iIndexGenome1+1,len(tListOfGenomeDbId)):
				sGenomeDbId2=tListOfGenomeDbId[iIndexGenome2]
				tListOfCouple.append((str(sGenomeDbId1),str(sGenomeDbId2)))
				
		chuckSize=len(tListOfCouple)/iMaxJobForParallelizeHLink
		iLasPacket=len(tListOfCouple)%iMaxJobForParallelizeHLink
		iIndex=0
		printDebug("----------Step 8 --3-- : Ending",sStdOutFile,bPrint)
		while iIndex<=iMaxJobForParallelizeHLink:
			tListOfGenomeId=[]
			for dbCouple in tListOfCouple[iIndex*chuckSize:(iIndex+1)*chuckSize]:
				tListOfGenomeId+=[str(dbCouple[0]),str(dbCouple[1])]
			print list(tListOfGenomeId) 
			if sMethodId is None:
				exit("sMethodId is None / 1")	
				sys.exit()
			#O2M:TODO:0.5
			ArgFile1.write(sVersionOrthocis+"\t"+str(bUSE_ONE_2_MANY)+"\t"+sMethodId+"\t"+conffileabs+"\t"+"\t".join(tListOfGenomeId)+"\n")
			iLineNumber1+=1
			iIndex+=1

		ArgFile1.close()
		printDebug("----------Step 8 --4-- : Ending",sStdOutFile,bPrint)
		ArgFile2=open(sTempName2,"w")

		#add fm 28 03
            
                for iIndexG1 in range(0,len(tListOfGenomeDbId)):
                        tListOfGenomeId=[]
                        tListOfGenomeId+=[str(tListOfGenomeDbId[iIndexG1])]
                       	if sMethodId is None:
                                exit("sMethodId is None / 2") 
				sys.exit()
			#O2M:TODO:0.5
			ArgFile2.write(sVersionOrthocis+"\t"+str(bUSE_ONE_2_MANY)+"\t"+sMethodId+"\t"+conffileabs+"\t"+"\t".join(tListOfGenomeId)+"\n")
			iLineNumber2+=1
		ArgFile2.close()		
		printDebug("----------Step 8 --5-- : Ending",sStdOutFile,bPrint)
		## activate this .? seei #TAG1
		doCleanMB=False
		if(doCleanMB):	
			sSQLinitHMember="delete hm from "+sVersionOrthocis+"_homology_member hm LEFT JOIN "+sVersionOrthocis+"_homology h on hm.homology_id=h.homology_id  where h.homology_id is NULL;"
			sSQLinitMember="delete m from "+sVersionOrthocis+"_member m LEFT JOIN "+sVersionOrthocis+"_homology_member hm on m.member_id=hm.member_id  where hm.member_id is NULL;"
			SubmitMySQLCommand(sSQLinitHMember,host,database,user,pw,db,bHaveMySQLdb,False)
			SubmitMySQLCommand(sSQLinitMember,host,database,user,pw,db,bHaveMySQLdb,False)
		sSQLinitGene="delete from "+sVersionOrthocis+"_gene ;alter table "+sVersionOrthocis+"_gene AUTO_INCREMENT=1;"
		print "delete gene"
		printDebug("----------Step 8 --6-- : Ending",sStdOutFile,bPrint)
		SubmitMySQLCommand(sSQLinitGene,host,database,user,pw,db,bHaveMySQLdb,False)
		sSQLinitHL="delete from "+sVersionOrthocis+"_homolink ;"
		print "delete homolink"
		SubmitMySQLCommand(sSQLinitHL,host,database,user,pw,db,bHaveMySQLdb,False)
		sSQLinitHL="delete from "+sVersionOrthocis+"_globalhomolinkstart ;"
                print "delete globalhomolinkstart"
                SubmitMySQLCommand(sSQLinitHL,host,database,user,pw,db,bHaveMySQLdb,False)

		#For all, analyze one gene against successively all other
		#Comparison are not oriented : no necessity to test g2 vs g1 if g1 vs g2 already made
		print tArrayOfArg
		#albane seconde zone commentee
	 	#nb genomes jobs: create tem tables
		printDebug("----------Homolink A",sStdOutFile,bPrint)
		printDebug("----------%s %s %s %s %s %s" % (sTempName2,iLineNumber2,jobexecdir+"/A",shomolinkA,sQueue,bHaveDrmaa),sStdOutFile,bPrint)
		TryParallelize(sTempName2,iLineNumber2,jobexecdir+"/A",shomolinkA,sQueue,bHaveDrmaa)
		#nb*nb genomes jobs: populate homilnk and genes

		#O2M:TODO:0.6
		printDebug("----------Homolink B",sStdOutFile,bPrint)
		TryParallelize(sTempName1,iLineNumber1,jobexecdir+"/B",shomolinkB,sQueue,bHaveDrmaa)
		#nb genomes jobs: delete temp tables	
		printDebug("----------Homolink C",sStdOutFile,bPrint)
		TryParallelize(sTempName2,iLineNumber2,jobexecdir+"/C",shomolinkC,sQueue,bHaveDrmaa)

		printDebug("----------Step 8 Fill Homolink table : Ending",sStdOutFile,bPrint)
		#DEBUG---------------------------------------------------
		DEBUG="""
		mysql -u [user] -p[pw] -h [host] [database] -e ""select * from [sVersionOrthocis]_gene limit 2;"
			+----------------+--------------------+--------------+
			| gene_id        | stable_id          | genome_db_id |
			+----------------+--------------------+--------------+
			| 12500001430054 | ENSPTRG00000019974 |          125 | 
			| 12500001434594 | ENSPTRG00000002715 |          125 | 
			+----------------+--------------------+--------------+
			
		mysql -u [user] -p[pw] -h [host] [database] -e ""select * from [sVersionOrthocis]_homolink limit 2;"
			+---------------+----------------+
			| first         | second         |
			+---------------+----------------+
			| 3300001426946 | 13500001426946 | 
			| 3300001427148 | 14900001427148 | 
			+---------------+----------------+
		"""
		#DEBUG---------------------------------------------------

		#O2M:

		ExecuteBashFile("python "+sdefineGlobalID+" "+conffile+" "+sVersionOrthocis);
	
	else:
		printDebug("----------Step 8 Fill Homolink table : No",sStdOutFile,bPrint)

	printDebug("FillHomolinkTable=false",sCurrentDebugFile,False)

else:
	printDebug("----------Step 3 Build Orthocis : No",sStdOutFile,bPrint)
	printDebug("----------Next Step : 9",sStdOutFile,bPrint)

printDebug("DoAllbase=false",sCurrentDebugFile,False)

###
if bRetrieveComponentConnex:
	cctablename=sVersionOrthocis+"_component_connex"
        globalcctablename=sVersionOrthocis+"_global_component_connex"
	printDebug("----------Step 9 Retrieve component connex : Starting",sStdOutFile,bPrint)
	if bCreateComponentConnexTable:
		printDebug("-Build component connex table : Starting",sStdOutFile,bPrint)		
		#O2M
        	#createComConnexTable(cctablename,"gene_id"," ismany=0 ",host,database,user,pw,db,sVersionOrthocis,bHaveMySQLdb)
        	if bUSE_ONE_2_MANY==True:
                	createComConnexTable(globalcctablename,"global_id"," ismany>=0 ",host,database,user,pw,db,sVersionOrthocis,bHaveMySQLdb)
		printDebug("-Build component connex table : Ending",sStdOutFile,bPrint)
	else:
		printDebug("-Build component connex table : No",sStdOutFile,bPrint)
		printDebug("CreateComponentConnexTable=false",sCurrentDebugFile,False)
	if bFillComponentConnexTable:
		printDebug("-Fill component connex table : Starting",sStdOutFile,bPrint)
		excludedSpecIDList=""
		if excludedGID and len(excludedGID)>1:
        		excludedSpecIDList+=""+excludedGID+"" # filter usefull when first steps already done will species including ortho errors
        	#executeCompoScript(sVersionOrthocis,host,database,user,pw,excludedSpecIDList,jobexecdir,scompco,sQueueCC,bHaveDrmaa,cctablename)
        	if bUSE_ONE_2_MANY==True:
			executeCompoScript(sVersionOrthocis,host,database,user,pw,excludedSpecIDList,jobexecdir,scompco,sQueueCC,bHaveDrmaa,globalcctablename,"GLOBAL")
		printDebug("-Fill component connex table : Ending",sStdOutFile,bPrint)
	else:
		printDebug("-Fill component connex table : No",sStdOutFile,bPrint)

	printDebug("FillComponentConnexTable=false",sCurrentDebugFile,False)

	printDebug("----------Step 9 Retrieve component connex : Ending",sStdOutFile,bPrint)
else:
	printDebug("----------Step 9 Retrieve component connex : No",sStdOutFile,bPrint)

printDebug("RetrieveComponentConnex=false",sCurrentDebugFile,False)

if bRetrieveClique:
	dico()
	printDebug("----------Step 10 Retrieve ortholog cliques : Starting",sStdOutFile,bPrint)
	if bCreateCliqueBoolTable:

		printDebug("-Build clique table : Starting",sStdOutFile,bPrint)
		createCliqueTable(sVersionOrthocis,host,database,user,pw,db,bHaveMySQLdb,dRelationGenomeDBId_SpecieName)
		if bUSE_ONE_2_MANY==True:
			createCliqueTable(sVersionOrthocis,host,database,user,pw,db,bHaveMySQLdb,dRelationGenomeDBId_SpecieName,"_global")

		printDebug("-Build clique table : Ending",sStdOutFile,bPrint)
	else:
		printDebug("-Build clique table : No",sStdOutFile,bPrint)

	printDebug("CreateCliqueBoolTable=false",sCurrentDebugFile,False)

	if bFillCliqueTables:
		if not bCreateCliqueBoolTable:
			printDebug("-Reset clique table",sStdOutFile,bPrint)
			sSQLreset2="DELETE FROM "+sVersionOrthocis+"_clique_element;"
			SubmitMySQLCommand(sSQLreset2,host,database,user,pw,db,bHaveMySQLdb,False)
			sSQLreset2="DELETE FROM "+sVersionOrthocis+"_global_clique_element;"
			SubmitMySQLCommand(sSQLreset2,host,database,user,pw,db,bHaveMySQLdb,False)

		printDebug("-Fill clique tables : Starting",sStdOutFile,bPrint)
		
		fillCliqueTable(sVersionOrthocis,host,database,user,pw,db,bHaveMySQLdb,iMaxJobForParallelize,conffileabs,jobexecdir,sDrmaaScript_RetrieveClique,sQueueClique,bHaveDrmaa)
		if bUSE_ONE_2_MANY==True:
			fillCliqueTable(sVersionOrthocis,host,database,user,pw,db,bHaveMySQLdb,iMaxJobForParallelize,conffileabs,jobexecdir,sDrmaaScript_RetrieveClique,sQueueClique,bHaveDrmaa,"_global")
		
		printDebug("-Fill clique tables : Ending",sStdOutFile,bPrint)
	else:
		printDebug("-Fill clique tables : No",sStdOutFile,bPrint)

	printDebug("FillCliqueTables=false",sCurrentDebugFile,False)

	printDebug("----------Step 10 Retrieve ortholog cliques : Ending",sStdOutFile,bPrint)
else:
	printDebug("----------Step 10 Retrieve ortholog cliques : No",sStdOutFile,bPrint)

printDebug("RetrieveClique=false",sCurrentDebugFile,False)

###
if bBuildHGNCtables:
	printDebug("----------Step 11 Build HGNC table : Starting",sStdOutFile,bPrint)
	
	if not os.path.isdir(sPathToHGNCfolder):
		exit("\nImpossible to find the folder : "+sPathToHGNCfolder+"\nProcess broken")
		sys.exit()
	else:
		tFiles=os.listdir(sPathToHGNCfolder)
		if "HGNC_complete_set.txt.gz" in tFiles:
			tCommand="gzip -d "+sPathToHGNCfolder+"/HGNC_complete_set.txt.gz"
			ExecuteBashFile(tCommand)
			sPathToHGNCfile=sPathToHGNCfolder+"/HGNC_complete_set.txt"
		elif "HGNC_complete_set.txt" in tFiles:
			sPathToHGNCfile=sPathToHGNCfolder+"/HGNC_complete_set.txt"
		else:
			exit("\nImpossible to find the file : "+sPathToHGNCfolder+"/HGNC_complete_set.txt (or txt.gz)\nProcess broken")
			sys.exit()
			
	printDebug("-Retrieve relation Human stableId/GeneId",sStdOutFile,bPrint)
	dRelationStableId_GeneId={}
	
	humangid=""
	#fm 03 2016
	sqlhsid="select genome_db_id from "+sVersionOrthocis+"_genome_db where name='homo_sapiens';"
	hsidtable=SubmitMySQLCommand(sqlhsid,host,database,user,pw,db,bHaveMySQLdb,True)
	for tRow in hsidtable:
                humangid=tRow[0]
	printDebug("  --"+"humangid:"+humangid,sStdOutFile,bPrint)
	#----
	sSQLinitHGNC="delete from "+sVersionOrthocis+"_asso_HGNC_clique ; delete from "+sVersionOrthocis+"_HGNC ;alter table "+sVersionOrthocis+"_HGNC AUTO_INCREMENT=1;"
	SubmitMySQLCommand(sSQLinitHGNC,host,database,user,pw,db,bHaveMySQLdb,False)

	sSQLretriveAllStableIdGeneIdInHuman="select stable_id,gene_id from "+sVersionOrthocis+"_gene where genome_db_id="+humangid+";"
	tSQLtable=SubmitMySQLCommand(sSQLretriveAllStableIdGeneIdInHuman,host,database,user,pw,db,bHaveMySQLdb,True)
	for tRow in tSQLtable:
		dRelationStableId_GeneId[tRow[0]]=tRow[1]
	iCurrentIndex=0		

	
	bFirst=True
 	printDebug("HGNCfile: "+str(sPathToHGNCfile),sStdOutFile,bPrint)	
	for sLine in open(sPathToHGNCfile):
		if bFirst:
			bFirst=False
			continue
		
		tLine=sLine.split("\t")
		sStableName=tLine[18]
		sOtherStableName=tLine[36]
		sFullName=tLine[2]
		sShortName=tLine[1]
		
		if sStableName=="":
			if sOtherStableName!="":
				sStableName=sOtherStableName
			else:
				continue
		
		try:
			sGeneId=dRelationStableId_GeneId[sStableName]
		except KeyError:
			print "KeyError:"+sStableName
			continue		
		
		iCurrentIndex+=1
		sSQLInsertInHGNC="insert into "+sVersionOrthocis+"_HGNC values ("+str(iCurrentIndex)+",'"+string2sql(sFullName)+"','"+string2sql(sShortName)+"');"
		SubmitMySQLCommand(sSQLInsertInHGNC,host,database,user,pw,db,bHaveMySQLdb,False)
		
	 
		sSQLretrieveCliqueId="select cm.clique_id from "+sVersionOrthocis+"_clique_element cm where cm.sp"+humangid+"="+sGeneId+";"
		tSQLtable=SubmitMySQLCommand(sSQLretrieveCliqueId,host,database,user,pw,db,bHaveMySQLdb,True)
		for tRow in tSQLtable:
			sCliqueId=tRow[0]
			sSQLInsertInHGNCassoClique="insert into "+sVersionOrthocis+"_asso_HGNC_clique values ("+str(iCurrentIndex)+","+sCliqueId+");"
			SubmitMySQLCommand(sSQLInsertInHGNCassoClique,host,database,user,pw,db,bHaveMySQLdb,False)
		
	printDebug("----------Step 11 Build HGNC table : Ending",sStdOutFile,bPrint)
else:
	printDebug("----------Step 11 Build HGNC table : No",sStdOutFile,bPrint)
	
printDebug("BuildHGNCtables=false",sCurrentDebugFile,False)

'''
if bCreateFasta:
	dico()
	bDoDeleteSeqTable=True
	printDebug("----------Step 12 Select Upstream in Fasta : Starting",sStdOutFile,bPrint)

	tSpeciesTable=SubmitMySQLCommand(sgenomeIDNSQLCmd,host,database,user,pw,db,bHaveMySQLdb,True)
	sTempName=getTempdir()+"TempArgFile"+str(random.random())+str(random.random())+".txt"
	ArgFile=open(sTempName,"w")

	iIndex=0
	while iIndex<len(tSpeciesTable):
		tLine=tSpeciesTable[iIndex]
		iIndex+=1
		spe=tLine[1]
		prfix="NA"
		if spe in species2data.keys() :
			if species2data[spe]['fasta_file_prefix'] is not None:
				prfix=species2data[spe]['fasta_file_prefix']
		ArgFile.write(sVersionOrthocis+"\t"+conffileabs+"\t"+tLine[0]+"\t"+spe+"\t"+fastaInDir+"\t"+fastaOutDir+"\t"+str(threshold)+"\t"+str(upstreamMode)+"\t"+biotype+"\t"+prfix+"\n")
		if bDoDeleteSeqTable:	
			sDelSeqTableSQL="delete from "+sVersionOrthocis+"_sequence_"+tLine[0]+";"
	     		SubmitMySQLCommand(sDelSeqTableSQL,host,database,user,pw,db,bHaveMySQLdb,False)
	
	ArgFile.close()
        #sTagCT=str(iIndex)+":"+str(upstreamMode)
	printDebug( "TryParallelize(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\") " %(sTempName,iIndex,jobexecdir+"/FastaBuild",sDrmaaScript_ExtractFromFasta,sQueueClique,bHaveDrmaa), sCurrentDebugFile,False)
	TryParallelize(sTempName,iIndex,jobexecdir+"/FastaBuild",sDrmaaScript_ExtractFromFasta,sQueueClique,bHaveDrmaa)

   			
	printDebug("----------Step 12 Select Upstream in Fasta : Ending",sStdOutFile,bPrint)
else:
	printDebug("----------Step 12 Select Upstream in Fasta : No",sStdOutFile,bPrint)

	printDebug("CreateFasta=false",sCurrentDebugFile,False)

printDebug("----------END----------",sStdOutFile,bPrint)
'''
##albane createFasta v2, on ne remplit que les tables fasta concernees par le speciesfilter
if bCreateFasta:

        bDoDeleteSeqTable=True
        printDebug("----------Step 12 Select Upstream in Fasta : Starting",sStdOutFile,bPrint)

        tSpeciesTable=SubmitMySQLCommand(sgenomeIDNSQLCmd,host,database,user,pw,db,bHaveMySQLdb,True)

        sTempName=getTempdir()+"TempArgFile"+str(random.random())+str(random.random())+".txt"
        ArgFile=open(sTempName,"w")

	filtered_species_list=speciesfilter.split(",")
        iIndex=0
	iFilteredSpeciesIndex=0
        while iIndex<len(tSpeciesTable):
                tLine=tSpeciesTable[iIndex]
                iIndex+=1
                spe=tLine[1]
		for species in filtered_species_list :
			p = re.compile(species)
			if p.match(spe) :
               			prfix="NA"
                		if spe in species2data.keys() :
                        		if species2data[spe]['fasta_file_prefix'] is not None:
                                		prfix=species2data[spe]['fasta_file_prefix']
                		ArgFile.write(sVersionOrthocis+"\t"+conffileabs+"\t"+tLine[0]+"\t"+spe+"\t"+fastaInDir+"\t"+fastaOutDir+"\t"+str(threshold)+"\t"+str(upstreamMode)+"\t"+biotype+"\t"+prfix+"\n")
				iFilteredSpeciesIndex+=1
                		if bDoDeleteSeqTable:   
                        		sDelSeqTableSQL="delete from "+sVersionOrthocis+"_sequence_"+tLine[0]+";"
                        		SubmitMySQLCommand(sDelSeqTableSQL,host,database,user,pw,db,bHaveMySQLdb,False)
        
        ArgFile.close()
        #sTagCT=str(iIndex)+":"+str(upstreamMode)
        printDebug( "TryParallelize(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\") " %(sTempName,iFilteredSpeciesIndex,jobexecdir+"/FastaBuild",sDrmaaScript_ExtractFromFasta,sQueueClique,bHaveDrmaa), sCurrentDebugFile,False)
        TryParallelize(sTempName,iFilteredSpeciesIndex,jobexecdir+"/FastaBuild",sDrmaaScript_ExtractFromFasta,sQueueClique,bHaveDrmaa)

                        
        printDebug("----------Step 12 Select Upstream in Fasta : Ending",sStdOutFile,bPrint)
else:
        printDebug("----------Step 12 Select Upstream in Fasta : No",sStdOutFile,bPrint)

        printDebug("CreateFasta=false",sCurrentDebugFile,False)

printDebug("----------END----------",sStdOutFile,bPrint)

