import sys
import shutil
import os
import re
import cPickle

from optparse import OptionParser


from ModuleFTP import *
from ModuleCreateTable_BuildOrthocis import *
from ModuleFonction_BuildOrthocis import *
from check_database import *

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
checkdir=oConfig.get('cf_checkdir')
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


check_generic_tables(host,database,user,pw,db,bHaveMySQLdb,checkdir)


	
