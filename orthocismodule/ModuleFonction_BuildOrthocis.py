import ConfigParser
import subprocess
import random
import os
import re
import glob
import time
import sys
import json


from ModuleNoDrmaa_BuildOrthocis import *

try:
	from ModuleDrmaa_BuildOrthocis import *
except:
	print "no drmaa available"

try:
        import MySQLdb
        locbHaveMySQLdb=True
except ImportError:
        locbHaveMySQLdb=False

modConf=None
useTempDir=True
tempDir=""

########################################################################
#CLASS
#Hack le fichier de conf. perl pour python
class FakeSecHead(object):
   def __init__(self, fp, sec):
      #section="mainconf"
      self.fp = fp
      self.sechead = '['+sec+']\n'
   def readline(self):
      if self.sechead:
         try: return self.sechead
         finally: self.sechead = None
      else: return self.fp.readline()

class CustomConfig(object):
	def __init__(self, confFile):
		self.file= confFile
		self.section= "mainconf"
		self.config= ConfigParser.ConfigParser()
		self.config.readfp(FakeSecHead(open(self.file),self.section))
	def get(self,key):
		return self.config.get(self.section,key)		

########################################################################
#FUNCTION
########################################################################

########################################################################



def trySubmitMySQLCommand(cmdSql,host,database,user,pw,db,bModule=False,bReturn=True, ctcount=0):
        ctmax=8
	ret=None
	sleeptime=10
        try:
           ret=SubmitMySQLCommand(cmdSql,host,database,user,pw,db,bModule=False,bReturn=True)
        except:
                print "warning error:", sys.exc_info()[0]
                if ctcount>ctmax:
                        raise
                else:
                        ctcount=ctcount+1
			time.sleep(sleeptime)
                        ret= trySubmitMySQLCommand(cmdSql,host,database,user,pw,db,bModule,bReturn, ctcount)
                        sleeptime=sleeptime+10
	return ret



def readJsonstring(s):
        s = s.replace("'", "\"")
        d = json.loads(s)
        return d


def populateGeneralConf(conffile):
	print "populateGeneralConf"
	oConfig = CustomConfig(conffile)
	global modConf
	modConf=oConfig
	return oConfig


#Progamm FUNCTION
#Translate string into boolean
def TranslateInBool(sString):
	sTemp=sString.lower()
	if sTemp=="true":
		return True
	if sTemp=="false":
		return False
	return TypeError	

def GetAllSequenceFromFasta(sFastaPath):
	dDicSeq={}
	sSeq=""
	sNameSeq="Empty"
	try:
		for sLine in open(sFastaPath):
			if sLine[0]==">":
				if sNameSeq!="Empty":
					dDicTempFasta[sNameSeq]=sSeq
					sSeq=""
				sNameSeq=sLine[1:].replace("\n","")
			else:
				sSeq+=sLine.replace("\n","")
				
		return(dDicSeq)
	except IOError:
		return IOError

#Print debug information
def printDebug(sString,sFile,bBool):
	FILE=open(sFile,"a")
	FILE.write(sString+"\n")
	FILE.close()
	if bBool:
		print(sString)


#Try Drmaa
def TryParallelize(pathArgFile,lineNumber,oudir,scriptName,queue,bBoolDrmaa,lg="python"):
	print "TryParallelize::1"
	if bBoolDrmaa:
		#LaunchDrmmaOn(pathArgFile,lineNumber,oudir,scriptName,queue)
		sys.exit("Don't use Drmaa for the moment...")
		print "TryParallelize::2"
	else:
		if len(oudir)>3:
			print "TryParallelize::3"+str(oudir)
			ExecuteBashCommand(["rm","-rf",oudir+"/*"])
			ExecuteBashCommand(["mkdir","-p",oudir])
		else:
			print "TryParallelize::4"
			sys.exit("Warning, path to short. Abort command rm -r on "+oudir)
		LaunchMultiQsubOn(pathArgFile,lineNumber,oudir,scriptName,queue,lg)
########################################################################
#MySQL FUNCTION
#Get the current id +1 in the table to simulate the auto-incrementation
def GetMaxGeneIdPlus1(sTable,sColumn,host,database,user,pw,db,bCheckMySQLdb):
	cmqSql="select max("+sColumn+") from "+sTable+";"
	iResult=0

	tSQLtable=SubmitMySQLCommand(cmqSql,host,database,user,pw,db,bCheckMySQLdb,True)
	for tRow in tSQLtable:
		if tRow[0]=="NULL":
			iResult=0
		else:
			iResult=int(tRow[0])
	
	iResult+=1
	return iResult

#Gestion of the MySQL command execution
def SubmitMySQLCommand(cmdSql,host,database,user,pw,db,bModule=False,bReturn=True):
	oResult=[]
        dosqllog=1
	print "-----------SubmitMySQLCommand------------"+cmdSql
	#if dosqllog==1:
		#print "-----------------------"
		#print cmdSql
		#print "-----------------------"	
	#
	if bModule:
		doConnectClose=False
		#print 1
		if db == None:
			#print "connect"
			doConnectClose=True
		        db = MySQLdb.connect(host=host,user=user,passwd=pw,db=database)
		#print 2
		cur = db.cursor()
		#print 3
		cur.execute(cmdSql)
		#print 4
		if bReturn:
			#print "fetch"
			for row in cur.fetchall():
				#print "app"
				oResult.append(row)
		#print 5
		cur.close()
		#print 6
		if doConnectClose:
			db.close()
		#print 7
	else:
		tArray=ExecuteMySQLCommand(cmdSql,host,database,user,pw)
		for row in tArray:
			oResult.append(row)
	return oResult


#Run a list of MysqlImport 
#WARNING ! File source name format : [X]_tableName.txt.gz

def ExecuteMultiMySQLImportFromAFolder(host,database,user,pw,path,bdDVersion,ensemblVersion,o2m,sspecGId=None):    
    #o2m: keep one 2 many orthologs 08/2017
	tFileList=os.listdir(path)
	for sFileName in tFileList:
		print("ExecuteMultiMySQLImportFromAFolder---"+sFileName+" ")
		doprocessfile=True

		 
		#if("gene.txt" in sFileName or "transcript.txt" in sFileName ):
		#	doprocessfile=True
		#	print "***************"+sFileName
		if doprocessfile and ".txt.gz" in sFileName and (not "_genome_db" in sFileName):

			
			sNewName=sFileName.replace(".gz","")
                        sTempName=sNewName.replace(".txt","")
			cutcmd=""
			dokeep=False
			if("gene_member" in sTempName):
				cutcmd="| cut -f1,2,6"
				dokeep=True
				sTempName=sTempName.replace("gene_member","member")
			if("homology_member" in sTempName):
                                cutcmd="| cut -f1,2"
				dokeep=True
			elif("homology" in sTempName):

				cutcmd="| cut -f1,2,3 |  grep "
				if o2m==True:
				   cutcmd+="ortholog_ "
				else:
				   cutcmd+="ortholog_one2one "

                                dokeep=True

#				ENUM:
#ens v84
#    ortholog_one2one
#    ortholog_one2many
#    ortholog_many2many
#    within_species_paralog
#    other_paralog
#    gene_split
#    between_species_paralog
#    alt_allele
#    homoeolog_one2one
#    homoeolog_one2many
#    homoeolog_many2many


				#doHomologyFull=False
				#if(doHomologyFull):
				#	sNewNamefull=sTempName+"_full.txt"
				#	fullcmd="zcat  "+path+"/"+sFileName+" > "+path+"/"+sNewNamefull
				#	ExecuteBashFile(fullcmd)
				#	sTableNameInBdD=bdDVersion+"_homology_full"
				#	ExecuteMySQLImportCommand(host,database,user,pw,path+"/"+sNewNamefull,sTableNameInBdD)		


			tCommand="zcat  "+path+"/"+sFileName+cutcmd+" > "+path+"/"+sNewName
			print ""+tCommand+"\n" 
			ExecuteBashFile(tCommand)
			sTableName=sTempName
			if ensemblVersion+"_" in sTableName:
				sTableName=sTableName.replace(ensemblVersion+"_","")
 	 
			sTableNameInBdD=bdDVersion+"_"+sTableName
			if sspecGId!=None:
				sTableNameInBdD=sTableNameInBdD+"_"+sspecGId
				print("-sspecGId!=None:"+str(sspecGId))
			else:
				print("WARNING!!-sspecGId==None!!")
			print("-"+sNewName+" import in database")
			ExecuteMySQLImportCommand(host,database,user,pw,path+"/"+sNewName,sTableNameInBdD)
			if not dokeep:
				ExecuteBashFile("rm -f "+path+"/"+sNewName)
		else:
			print sFileName+ " not processed"
#SQL command to the special mysql load data	syntax
def ExecuteMySQLImportCommand(host,database,user,pw,path,table):
        cmdSql= "LOAD DATA LOCAL INFILE  \\\""+path+"\\\" REPLACE INTO TABLE  "+table+";"  
	mysql=["mysql", " --max_allowed_packet=1G ","--batch", "-u",user+"","-p"+pw,"-h",host+"",database+"","-e","\""+cmdSql+"\""]
	myssqlTemp=" ".join(mysql)
	print "ExecuteMySQLImportCommand:cmd:"
	print "-------------------"
	print myssqlTemp
	print "-------------------"
 	return ExecuteBashFile(myssqlTemp)

#Make a SQL command
def ExecuteMySQLCommand(cmdSql,host,database,user,pw, outfile=[]):
	cmdSql=cmdSql.replace("\n"," ")
	mysql=["mysql"," --max_allowed_packet=1G ", "--batch", "-u",user+"","-p"+pw,"-h",host+"",database+"","-e","\""+cmdSql+"\""]+outfile
	myssqlTemp=" ".join(mysql)	
	outp=ExecuteBashFile(myssqlTemp,True,'Lost\sconnection\sto\sMySQL\sserver')
 	return outp 


def ExecuteMySQLCommandFromFile(cmdSql,host,database,user,pw):

	cmdSql=cmdSql.replace("\n"," ")

	cmdSqlFile=getTempdir()+"Tempsql"+str(random.random())+str(random.random())+".sql"

	TMPF=open(cmdSqlFile,"w")
	TMPF.write(cmdSql)
	TMPF.close()


	mysql=["mysql", "--batch", "-u",user+"","-p"+pw,"-h",host+"",database+"","< ","\""+cmdSqlFile+"\""]
	myssqlTemp=" ".join(mysql)	

 	ret= ExecuteBashFile(myssqlTemp)
	os.remove(cmdSqlFile)
	return ret	

########################################################################
#FILE MANIPULATION FUNCTION
#Two non-exclusive possibility for the fasta : primary_assembly and top-level. If both, priority to assembly
def generatePrimOrTopFastaFromArchive(inpath,version,outpath):
	dIsPresentPrimary=False
	dIsPresentTopLevel=False
	iArePresentBoth=0
	tFileList=os.listdir(inpath)

	spec2n={}
	for sArch in tFileList:
		if "toplevel" in sArch and ".gz" in sArch:
			specn=sArch.split(".")[0]
			spec2n[specn]=sArch
	for sArch in tFileList:
		if "primary" in sArch and ".gz" in sArch:
			specn=sArch.split(".")[0]
			spec2n[specn]=sArch

 
 	for specname in spec2n:
		soFileName=spec2n[specname]
		sNewName=soFileName.replace(".fa.gz","")
		tCommand="zcat "+inpath+"/"+soFileName +"> "+outpath+"/"+sNewName+".fasta"
		print tCommand
		ExecuteBashFile(tCommand)
 
 

#Rename file according to the bdd convention
def RenameAllFileTxt(path,version):
	tFileList=os.listdir(path)
	for sFile in tFileList:
		if ".txt.gz" in sFile:
			os.renames(path+"/"+sFile,path+"/"+version+"_"+sFile)
	return True

#Gunzip the file genome_db in ensembl_compara folder
def GetCommandGunzip_genomeDB(path,toString_version):
	sFilePath=path+"/ensembl_compara_"+toString_version
	
	tCommand="zcat "+sFilePath+"/*genome_db.txt.gz  > "+sFilePath+"/"+"genome_db.txt "
 

	ExecuteBashFile(tCommand)
#################

def getTempdir():
	global tempDir
	global useTempDir
	global modConf
	dpath=""
	if useTempDir:
		if(tempDir==None or tempDir==""):
			tempDir= modConf.get('cf_tempdir')
		dpath=tempDir+"/"
	return dpath
########################################################################
#BASH EXECUTION FUNCTION
def ExecuteBashFile(scriptfile,bFirstLine=True, msgretry="", rct=0):
	sName=getTempdir()+"TempExceBF"+str(random.random())+str(random.random())+".sh"
	FileTemp=open(sName,"w")
	FileTemp.write(scriptfile)
	FileTemp.close()

	tTable=ExecuteBashCommand(["bash",sName])

	os.remove(sName)

	if tTable[0]==0:
		tRows=[]
		sOut=tTable[1]
		tOut=sOut.split("\n")
		for sLine in tOut:
			if not bFirstLine:
				if len(sLine)!=0:
					sLine=sLine.strip()
					tRows.append(re.split("\t",sLine))
			else:
				bFirstLine=False
		return tRows
	else:
		errmsg=tTable[2]
		if  msgretry!="" and rct < 4 :
			p = re.compile(".*"+msgretry)
			m = p.match(errmsg)
			if m:
				rct=rct+1
				print errmsg+ "\nfailover policy, going to retry  ("+str(rct)+") "
				time.sleep(60)			
                        	return ExecuteBashFile(scriptfile,bFirstLine, msgretry, rct)
		
		print errmsg
                sys.exit(scriptfile+" exec command error: "+errmsg  )			


#Make a Bash command
def ExecuteBashCommand(cmdArray):
	sp = subprocess.Popen(cmdArray, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = sp.communicate()
	if err:
	    print "standard error of subprocess:"
	    print err
	    if out:
            	print "standard output of subprocess:"
            	print out	
	print "returncode of subprocess:", sp.returncode
	return [sp.returncode,out,err]

########################################################################
#CREATE CONTENT FILE
def TemplateRunDrmaaScript(pythonScriptName,ArrayOfArrayofArg,workDir,outDir,queue):

	localdir=os.getcwd()
	modulePath="import sys\nsys.path.append('"+localdir+"')\n"

	#Stock Array of Arg
	bFirst=True
	sFileArrayPath=workDir+"/ArrayOfArg.txt"
	FileTemp=open(sFileArrayPath,"w")
	for tArray in ArrayOfArrayofArg:
		if bFirst:
			bFirst=False
		else:
			FileTemp.write("\n")
		FileTemp.write("\t".join(str(item) for item in tArray))
	FileTemp.close()

	return modulePath+"from ModuleDrmaa_BuildOrthocis import *\nLaunchDrmaa('"+pythonScriptName+"','"+sFileArrayPath+"','"+workDir+"','"+outDir+"','"+queue+"')"


########################################################################
#DRMAA LAUNCHER FUNCTION
def LaunchDrmmaOn(ArrayOfArg,oudir,scriptName,queue):
	drmaaFolder=oudir+"/drmaaTempFolder"
	try:
		os.makedirs(drmaaFolder)
	except OSError :
		print(drmaaFolder+" already exist")

	sRunDrmaaScripContent=TemplateRunDrmaaScript(scriptName,ArrayOfArg,drmaaFolder,drmaaFolder,queue)
	sNameRunDrmaaScript="Temp_RunDrmaa.py"

	FileTemp=open(drmaaFolder+"/"+sNameRunDrmaaScript,"w")
	FileTemp.write(sRunDrmaaScripContent)
	FileTemp.close()

	sBashScriptContent="""
	export DRMAA_LIBRARY_PATH=/usr/local/sge/lib/lx24-amd64/libdrmaa.so
	source /local/env/envpython-2.6.4.sh
	python """+drmaaFolder+"/"+sNameRunDrmaaScript
	ExecuteBashFile(sBashScriptContent)
	
	return True



def string2sql(s):
	col  = str(s).replace("'", "\\'")
	return col


def insensitive_glob(pattern):
    def either(c):
        return '[%s%s]'%(c.lower(),c.upper()) if c.isalpha() else c
    return glob.glob(''.join(map(either,pattern)))

def filterSpecFastaFromDir(indir,species,sfx,prefix="NA"):
	if prefix=="NA":
		listing = insensitive_glob(indir+"*"+species+"*"+sfx)
	else:
		listing = insensitive_glob(indir+prefix+"*"+sfx)	
	tab=[]	
	for filename in listing:
		tab.append(filename)
	return tab

def executeCompoScript(sVersionOrthocis,host,database,user,pw,excludedSpecIDList,jobexecdir,scompco,sQueueCC,bHaveDrmaa,cctablename,globalPrefix=""):
	hltabfile=getTempdir()+globalPrefix+"HL.tab"
	sqlglt=""
        if globalPrefix != "":
		sqlglt="_global"
	ExecuteMySQLCommand("select * from "+sVersionOrthocis+sqlglt+"_homolink;",host,database,user,pw, ["|","sed","'1d'",">",hltabfile])
	cctabfile=getTempdir()+globalPrefix+"CC.cc"
	sTempNameCC=getTempdir()+globalPrefix+"ccarg.txt"
	ArgFile=open(sTempNameCC,"w")
	ArgFile.write(hltabfile+"\t"+globalPrefix+"cc_count.txt\t"+excludedSpecIDList+"\t>\t"+cctabfile)
	ArgFile.close()
	#generate components connex
	TryParallelize(sTempNameCC,1,jobexecdir+"/CC"+globalPrefix,scompco,sQueueCC,bHaveDrmaa,"perl")
	#import cc in db
	ExecuteMySQLImportCommand(host,database,user,pw,cctabfile,cctablename)

	
def createCliqueTable(sVersionOrthocis,host,database,user,pw,db,bHaveMySQLdb,dRelationGenomeDBId_SpecieName,globalTag=""):
	sSQLcreateTable="drop table if exists "+sVersionOrthocis+globalTag+"_clique_element;"
	sSQLcreateTable+="create table "+sVersionOrthocis+globalTag+"_clique_element (clique_id bigint(14)"
	sSQLcreateTable+=" AUTO_INCREMENT, primary key (clique_id));"
	SubmitMySQLCommand(sSQLcreateTable,host,database,user,pw,db,bHaveMySQLdb,False)
	
	tAllGenomeDbId=sorted(dRelationGenomeDBId_SpecieName.values())
	for iGenomeDbId in tAllGenomeDbId:
		sSQLaddColumn="alter table "+sVersionOrthocis+globalTag+"_clique_element add sp"+str(iGenomeDbId)+" bigint(14) default NULL ;"
		SubmitMySQLCommand(sSQLaddColumn,host,database,user,pw,db,bHaveMySQLdb,False)

def fillCliqueTable(sVersionOrthocis,host,database,user,pw,db,bHaveMySQLdb,iMaxJobForParallelize,conffileabs,jobexecdir,sDrmaaScript_RetrieveClique,sQueueClique,bHaveDrmaa,globalTag=""):
	sSQLretrieveAllComponentId="select distinct component_id from "+sVersionOrthocis+globalTag+"_component_connex where component_id is not null;"

	sTempName=getTempdir()+"TempArgFile"+str(random.random())+str(random.random())+".txt"
	ArgFile=open(sTempName,"w")
	iLineNumber=0

	tSQLtable=SubmitMySQLCommand(sSQLretrieveAllComponentId,host,database,user,pw,db,bHaveMySQLdb,True)
	 
	iPacket=len(tSQLtable)/iMaxJobForParallelize
	iLasPacket=len(tSQLtable)%iMaxJobForParallelize

	iIndex=0
	while iIndex<=iMaxJobForParallelize:
		tListOfGeneId=[]
		for tUniqGeneId in tSQLtable[iIndex*iPacket:(iIndex+1)*iPacket]:
			tListOfGeneId+=tUniqGeneId
	 
		ArgFile.write(sVersionOrthocis+"\t"+conffileabs+"\t"+globalTag+"\t"+"\t".join(tListOfGeneId)+"\n")
		iIndex+=1
		iLineNumber+=1
	ArgFile.close()
	 
	TryParallelize(sTempName,iLineNumber,jobexecdir+"/CL",sDrmaaScript_RetrieveClique,sQueueClique,bHaveDrmaa)





