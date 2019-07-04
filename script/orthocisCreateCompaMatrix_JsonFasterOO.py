import random
from optparse import OptionParser
import re
import ConfigParser
import subprocess
import time
import json
import math
from xml.dom.minidom import parseString
import ast
from OrthocisHomology import *
from pprint import pprint 
import os

sVersionScript="v30"
##################
iHomogeneityLimit=750000
# iHomogeneityLimit is a tehcniocal fix for performances
#justification:
# we considere than cis regions of genes with a very large nulbre of hits are junk
###################
#HOMOLOGY_SCRIPT_MODE="python"
HOMOLOGY_SCRIPT_MODE="cpp"
########################################################################

class FinalDataStruct:
    """Class 
         related to final results
    """
    
    def __init__(self, arg): 

        self.sHost  = arg.sHost
	self.sDatabase = arg.sDatabase
	self.sUser = arg.sUser
	self.sPw = arg.sPw
        self.sVersionDb = arg.sVersionDb
        self.dSpeciesCount = dict()
        self.internal_genes = list()        
        self.genes = list()

    def computeGenes(self,limit, arg):
	if len(self.internal_genes) >0:
		sSQLcommand="select stable_id,genome_db_id from "+self.sVersionDb+"_gene"
           	sSQLcommand+=" where gene_id in ("+",".join(self.internal_genes)+") limit "+str(limit)+";"
           	#print sSQLcommand
           	#tSQLtable=SubmitMySQLCommand(sSQLcommand,sHost,sDatabase,sUser,sPw)
		query = SQLquery(sSQLcommand, arg)
           	for tRows in query.result :
			sENSid=tRows[0]
			self.genes.append(sENSid)

def exportDataDictHom(dDico,sTempFileName="Output.txt"):
	FILE=open(sTempFileName,"w")
	FILE.write("Gene\tHit\tStart\tEnd\tPvalue\tStrand\tSequence\n")
	for sGene in dDico:
		for sHit in dDico[sGene]:
			FILE.write(sGene+"\t"+sHit+"\t"+dDico[sGene][sHit]['Start']+"\t"+dDico[sGene][sHit]['End']+"\t"+dDico[sGene][sHit]['Pvalue']+"\t"+dDico[sGene][sHit]['Strand']+"\t"+dDico[sGene][sHit]['Sequence']+"\n")
	FILE.close()

def exportDataList(tExceptedPos,sTempFileName="Output2.txt"):
	FILE=open(sTempFileName,"w")
 
	for el in tExceptedPos:
		FILE.write(el+"\n")
	FILE.close()


def ExecuteBashCommand(cmdArray):
        #print(str(cmdArray))
	sp = subprocess.Popen(cmdArray, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = sp.communicate()
	if err:
	    print "standard error of subprocess:"
	    print err
	    if out:
            	print "standard output of subprocess:"
            	print out	
	#print "returncode of subprocess:", sp.returncode
	return [sp.returncode,out,err]

 
def launchHomology4MultiHits(mode,dWorkGene2Data,tExceptedPosition,fIdentityMax):
	if mode=="python":
		 return homology4MultiHits(dWorkGene2Data,tExceptedPosition,fIdentityMax)
	elif mode=="cpp":
		 
		try:
			res='''
			#example output
			__RESULT__
			<result>
			<score>
			31.0
			</score>
			<hits>
			{'15800000004639': {'49057241': {'Start': '-177440', 'End': '-177425', 'Pvalue': '8.7e-06', 'Strand': '1', 'Sequence': 'TGACCTTAAGTGATCC'}}, '15800000004640': {'49057224': {'Start': '-177540', 'End': '-177425', 'Pvalue': '8.1e-06', 'Strand': '1', 'Sequence': 'TTACCTTAAGTGATCC'}}}
			</hits>
			</result>

			'''

			sName1="datatmp1"+str(random.random())+".txt"
			sName2="datatmp2"+str(random.random())+".txt"
			exportDataDictHom(dWorkGene2Data,sName1)
			exportDataList(tExceptedPosition,sName2)
		        binloc=os.environ["HOMOGEN"]	
			if binloc == None or len(binloc)<1:
				sys.exit("Error : env var HOMOGEN not defined")	 
                        cmd=["bash","-c", ""+binloc+"/"+"OrthocisHomology"+" "+sName1+" "+sName2+" "+str(int(fIdentityMax)) +"" ]
                        #pprint(cmd)
			tTable=ExecuteBashCommand(cmd)
			os.remove(sName1)
			os.remove(sName2)
			#print tTable
 			rcode=tTable[0]
			rout=tTable[1]
			rerr=tTable[2]
			res=rout
			rlist=res.split('__RESULT__')
			resu=rlist[-1]
			xmldoc = parseString(resu) 
			itemlist = xmldoc.getElementsByTagName('result')
			scorenode=itemlist[0].getElementsByTagName('score')[0]
			hitsnode=itemlist[0].getElementsByTagName('hits')[0]
			fIdentityScore=float(scorenode.firstChild.wholeText.strip())
			dBestHit4Genes=str(hitsnode.firstChild.wholeText.strip())
			
			dBestHit4Genes= ast.literal_eval(str(hitsnode.firstChild.wholeText.strip() ))
			#print fIdentityScore
			#print dBestHit4Genes
			return [fIdentityScore,dBestHit4Genes]
		except Exception, e:
		  print str(e)	
		  return None


########################################################################
#Class
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

def populateGeneralConf(conffile):
	oConfig = CustomConfig(conffile)
	return oConfig

########################################################################
#Function

def StepIncrease(iStepValue,sString,iDepth):
	sLine=""+"\t"*iDepth
	sLine+="- Step "+str(iStepValue)
	sLine+=" : "+sString
	print sLine
	return (iStepValue+1)

def GetDifferenceBetweenList(tRef,tComp):
	setComp=set(tComp)
	tResult=[element for element in tRef if element not in setComp]
	return tResult

def GetCommonValue(dDico,tExcepted,tTargetedSpecies):
	#Rebuild list of sequences
	tListOfSequences=[]
	for sGeneId in dDico:
		if sGeneId[:-11] in tTargetedSpecies:
			for sScanId in dDico[sGeneId]:
				sNewSeq=""
				for iIndex in range(len(dDico[sGeneId][sScanId]["Sequence"])):
					if iIndex not in tExcepted:
						sNewSeq+=dDico[sGeneId][sScanId]["Sequence"][iIndex]
				tListOfSequences.append(sNewSeq)

	#Compute common column
	tConsensus=list(tListOfSequences[0])
	for sSeq in tListOfSequences[1:]:
		for iIndex in range(len(tConsensus)):
			if sSeq[iIndex]!=tConsensus[iIndex]:
				tConsensus[iIndex]="*"
	
	return len([X for X in tConsensus if X!="*"])

#Gestion of the MySQL command execution
def SubmitMySQLCommand(cmdSql,host,database,user,pw):
	cmdSql = cmdSql.replace("\n"," ")
	mysql=["mysql", "--batch", "-u",user+"","-p"+pw,"-h",host+"",database+"","-e","\""+cmdSql+"\""]
	myssqlTemp=" ".join(mysql)	
 	return ExecuteBashFile(myssqlTemp)

def ExecuteBashFile(scriptfile,bFirstLine=True):
	sName="TempCompa"+str(random.random())+".sh"
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
		print tTable[2]
		exit("Command SQL false: "+scriptfile)


#Oriented object transformation

class Arguments : 
	def __init__(self) :

		parser = OptionParser(conflict_handler="resolve") #resolve : the default -h/--help option is disabled in favor of -h/--hostname
 
		parser.add_option("-c","--conf", dest="conf")
		parser.add_option("-h","--hostName", dest="hostName")
		parser.add_option("-u","--user", dest="user")
		parser.add_option("-b","--base", dest="base")
		parser.add_option("-p","--password", dest="password")

		parser.add_option("-a","--argFile", dest="argFile")
		parser.add_option("-j","--jsonFilePath", dest="jsonFilePath")
		parser.add_option("-s","--stringJsonData", dest="stringJsonData")

		(options, args) = parser.parse_args()
		
		self.parser = parser

		self.options = options
		
		self.sAbsentValue = "AbsentValue"

	def which_varsource(self, options) :
		options = self.options
		sArgFile = self.options.argFile
		bHaveArgFile = True
		bHaveJsonFile = False
		self.ArgFile = sArgFile
		self.source = self.ArgFile
		if not sArgFile:
			bHaveArgFile = False
			sJsonFilePath = self.options.jsonFilePath
			self.JsonPath = sJsonFilePath
			self.source = self.JsonPath
			if not sJsonFilePath:
				bHaveJsonFile=False
				sJsonStringData = self.options.stringJsonData
				self.JsonString = sJsonStringData
				self.source = self.JsonString
				if not sJsonStringData:
					sys.exit("Error : no Options -a or -j or -s defined, process broken")
	
		#check si on est en mode chaine, conffile ou path

	def getMatrixSize(self) :  #step0
		
		sSQLcommand="select motif_id from "+self.sVersionDb+"_job where job_id="+str(self.sJobId)+";"
		query = SQLquery(sSQLcommand,self)

		if len(query.result)  <= 0 :
			sSQLcommand="select motif_id,job_id from "+self.sVersionDb+"_job ;"
			query = SQLquery(sSQLcommand,self)
			print(str(query.result))
			exit("Error : no motif id for sJobId=%s" %(self.sJobId))

		self.sMatrixId = query.result[0][0]

		sSQLcommand="select size from "+self.sVersionDb+"_motif where motif_id="+self.sMatrixId+";"
		query = SQLquery(sSQLcommand,self)
		sMatrixSize=query.result[0][0]
		iMatrixSize=int(sMatrixSize)


	def add_arguments(self, source) :
		#selon la variable de check, lis et charge les arguments en consequence (etape de domain=1)
		self.sRegion=None
		if source == self.ArgFile :
			oArgs =CustomConfig(sArgFile)
			self.sGetResult=oArgs.get(sJsonKeyGetResult)
			sJobId=oArgs.get(sJsonKeyJobId);
		
			sRefSpecies=oArgs.get(sJsonKeyRefSpecies);
			sMatchSpecies=oArgs.get(sJsonKeyMatchSpecies);
			sOrthoSpecies=oArgs.get(sJsonKeyOrthoSpecies);
			sDisplaySpecies=oArgs.get(sJsonKeyDisplaySpecies);
		
			sScoreValueMin=oArgs.get(sJsonKeyValueMin);
			sScoreValueMax=oArgs.get(sJsonKeyValueMax);
	
			sIdentityFilter=oArgs.get(sJsonKeyHomologyFilter)
			#sMatrixSize=oArgs.get(sJsonKeyMatrixSize)
			sIdentityMin=oArgs.get(sJsonKeyHomologyMin)
			sIdentityMax=oArgs.get(sJsonKeyHomologyMax)
			sSingleHitAnalysis=oArgs.get(sJsonWorkOnHit)
			sExceptedPosition=oArgs.get(sJsonKeyHomologyExceptedPosition)
			self.sRegion=oArgs.get(sJsonKeyRegion)	

		elif source == self.JsonPath :
			dJsonData=json.load(open(sSpeciesJsonPath))
		else :
			sJsonStringData=self.JsonString.replace("'","\"")
			dJsonData=json.loads(sJsonStringData)
		print dJsonData
		self.sGetResult=str(dJsonData[self.sJsonKeyFilter][self.sJsonKeyGetResult])
		self.sJobId=str(dJsonData[self.sJsonKeyFilter][self.sJsonKeyJobId])
	
		self.sRefSpecies=str(dJsonData[self.sJsonKeySpeciesId][self.sJsonKeyRefSpecies])
		self.sMatchSpecies=str(dJsonData[self.sJsonKeySpeciesId][self.sJsonKeyMatchSpecies])
		self.sOrthoSpecies=str(dJsonData[self.sJsonKeySpeciesId][self.sJsonKeyOrthoSpecies])
		self.sDisplaySpecies=str(dJsonData[self.sJsonKeySpeciesId][self.sJsonKeyDisplaySpecies])
	
		self.sSingleHitAnalysis=str(dJsonData[self.sJsonKeyFilter][self.sJsonToolConf][self.sJsonWorkOnHit])
		self.sScoreValueMin=str(dJsonData[self.sJsonKeyFilter][self.sJsonToolConf][self.sJsonKeyValueMin])
		#self.otherGenes = otherGenes
		#self.otherGenes = otherGenes
		self.sScoreValueMax=str(dJsonData[self.sJsonKeyFilter][self.sJsonToolConf][self.sJsonKeyValueMax])
		#Hack-> Json = Web interface = homology ever checked during script processing
		self.sIdentityFilter=str(dJsonData[self.sJsonKeyFilter][self.sJsonToolConf][self.sJsonKeyHomologyFilter])
		
		#print("@@"+str(dJsonData[sJsonKeyFilter][sJsonToolConf][sJsonKeyHomologyFilter]))   
		#sIdentityFilter="True"
		#sMatrixSize=str(dJsonData[sJsonToolConf][sJsonKeyMatrixSize])
		sIdentityMin=str(dJsonData[self.sJsonKeyFilter][self.sJsonToolConf][self.sJsonKeyHomologyMin])
		sIdentityMax=str(dJsonData[self.sJsonKeyFilter][self.sJsonToolConf][self.sJsonKeyHomologyMax])
		sExceptedPosition=str(dJsonData[self.sJsonKeyFilter][self.sJsonToolConf][self.sJsonKeyHomologyExceptedPosition])
		self.sRegion=dJsonData[self.sJsonKeyFilter][self.sJsonToolConf][self.sJsonKeyRegion]
	
		if self.sRegion is not None :
	
			if self.sRegion=='upstream':
				self.sregionClause=" and region=1 "
			if self.sRegion=='genefrag':
				self.sregionClause=" and region=2 "
			if self.sRegion=='upstream+genefrag':
				self.sregionClause=" and region in (1,2)"

	def add_database_parameters(self) :
	
		#selon la variable de check, lis et charge les arguments en consequence (etape de domain=1)
		if self.source == self.ArgFile :
			confile = self.options.conf
			oConfig = populateGeneralConf(conffile)
			self.sHost  = oConfig.get('host');
			self.sDatabase = oConfig.get('database');
			self.sUser = oConfig.get('user');
			self.sPw = oConfig.get('pw');
		else :
	
			self.sHost  = self.options.hostName
			self.sDatabase = self.options.base
			self.sUser = self.options.user
			self.sPw = self.options.password
			if not self.sHost:
				sys.exit("Error : no database login hostName -h  defined, process broken")
			if not self.sDatabase:
				sys.exit("Error : no database login baseName -b  defined, process broken")
			if not self.sUser:
				sys.exit("Error : no database login user -u  defined, process broken")
			if not self.sPw:
				sys.exit("Error : no database login password -p  defined, process broken")

	def check_validity(self) :
		#verifie si les options sont valides
		try:
			int(self.sJobId)
		except:
			exit("Error : No JobId defined, process broken")

		print self.sGetResult

		if self.sGetResult.lower()=="true" or self.sGetResult=="1":
			self.bGetResult=True
			print "do get result"
		elif self.sGetResult.lower()=="false" or self.sGetResult=="0" or self.sGetResult==sAbsentValue:
			self.bGetResult=False
			print "do count"
		else:
			exit("Error : GetResult option must be a boolean, process broken")

		if self.sRefSpecies=="":
			exit("Error : No RefSpecies defined, process broken")
		try:
			int(self.sRefSpecies)
		except ValueError:
			exit("Error : RefSpecies must be defined by integer, process broken")
		if self.sMatchSpecies=="":
			self.tMatchSpecies=[]
			print("Warning : no MatchSpecies defined")
		else:
			try:
				self.tMatchSpecies=self.sMatchSpecies.split(self.sOptionInternalItemSeparator)
				[int(X) for X in self.tMatchSpecies]
			except ValueError:
				exit("Error : MatchSpecies must be defined by integer, process broken")
		if self.sOrthoSpecies=="":
			self.tOrthoSpecies=[]
			print("Warning : no OrthoSpecies defined")
		else:
			try:
				self.tOrthoSpecies=self.sOrthoSpecies.split(self.sOptionInternalItemSeparator)
				[int(X) for X in self.tOrthoSpecies]
			except ValueError:
				exit("Error : OrthoSpecies must be defined by integer, process broken")
		if self.sDisplaySpecies=="":
			self.tDisplaySpecies=[]
			print("Warning : no DisplaySpecies defined")
		else:
			try:
				self.tDisplaySpecies=self.sDisplaySpecies.split(self.sOptionInternalItemSeparator)
				[int(X) for X in self.tDisplaySpecies]
			except ValueError:
				exit("Error : DisplaySpecies must be defined by integer, process broken")


	def add_constants(self) :
		#ajoute les constantes type item_separator, etc
		#CONSTANT
		self.sVersionDb="current"
		self.sJsonOutputName="countOO.json"
		self.sOptionInternalItemSeparator="-"

		self.sJsonKeySpeciesId="SpeciesId"
		self.sJsonKeyRefSpecies="RefSpecies"
		self.sJsonKeyMatchSpecies="MatchSpecies"
		self.sJsonKeyOrthoSpecies="OrthoSpecies"
		self.sJsonKeyDisplaySpecies="DisplaySpecies"

		self.sJsonKeyFilter="Filter"
		self.sJsonKeyGetResult="GetResult"
		self.sJsonKeyJobId="JobId"

		self.sJsonToolConf="toolconf"
		self.sJsonWorkOnHit="WorkOnAllHit"
		#sJsonKeyScoreFilter="ScoreFilter"
		self.sJsonKeyValueMin="ValueMin"
		self.sJsonKeyValueMax="ValueMax"
		self.sJsonKeyHomologyFilter="HomologyFilter"
		self.sJsonKeyHomologyOn="HomologyOn"
		#sJsonKeyMatrixSize="MatrixSize"
		self.sJsonKeyHomologyMin="HomologyMin"
		self.sJsonKeyHomologyMax="HomologyMax"
		self.sJsonKeyHomologyExceptedPosition="exceptedPosition"
		self.sJsonKeyRegion="scanRegion"

	def concatenate(self) :
		#concatene toutes les options dans une variable
		#Concatenate All species in one variable
		self.tAllSpecies=[self.sRefSpecies]+self.tMatchSpecies+self.tOrthoSpecies+self.tDisplaySpecies
		self.tCheckedSpecies=[self.sRefSpecies]+self.tMatchSpecies

		self.bApplyScoreFilter=False
		bApplyMinScore=False
		bApplyMaxScore=False
		if self.sScoreValueMin=="":
			sScoreValueMin="0"
			print("Warning : no ScoreValueMin defined, default value 0 used")
		else:
			try:
				self.fMinScore=float(self.sScoreValueMin)
				self.bApplyMinScore=True
				self.bApplyScoreFilter=True
			except ValueError:
				exit("Error : ValueMin must be defined by float, process broken")
		if self.sScoreValueMax=="":
			self.sScoreValueMax="Nothing"
			print("Warning : no ScoreValueMax defined, no value used")
		else:
			try:
				self.fMaxScore=float(self.sScoreValueMax)
				self.bApplyMaxScore=True
				self.bApplyScoreFilter=True
			except ValueError:
				exit("Error : ValueMax must be defined by float, process broken")

		self.bSingleHitAnalysis=False
		self.tExceptedPosition=[]
		if self.sSingleHitAnalysis=="all":
			print("WorkOnHit option : all (all avalaible hits on a gene are considered)")
			self.bSingleHitAnalysis=False
		elif self.sSingleHitAnalysis=="best":
			print("WorkOnHit option : best (only the best hit on a gene is considered)")
			self.bSingleHitAnalysis=True
		else:
			exit("Error : WorkOnHit must be defined by 'all' or 'best', process broken")

		self.bCheckHomogeneous=False
		if self.sIdentityFilter=="" or self.sIdentityFilter.lower()=="false":
			print("Warning : IdentityFilter fixed to False, no homogeneous checking will be made")
		elif len(self.tMatchSpecies+self.tOrthoSpecies+self.tDisplaySpecies)==0:
			print("Warning : only one species in input, no homogeneity checking will be made")
		else:
			self.bCheckHomogeneous=True
	
			#if sMatrixSize=="":
				#exit("Error : MatrixSize must be defined by integer, process broken")
			#try:
				#iMatrixSize=int(sMatrixSize)
			#except ValueError:
				#exit("Error : MatrixSize must be defined by integer, process broken")
	
			try:
				self.fIdentityMin=float(sIdentityMin)
			except ValueError:
				print("Warning : HomologyMin must be defined by float, default value 0 used")
				self.fIdentityMin=0
			try:
				self.fIdentityMax=float(sIdentityMax)
			except ValueError:
				exit("Error : HomologyMax must be defined by float, process broken")

			if self.sExceptedPosition!="" and self.sExceptedPosition!=sAbsentValue and self.sExceptedPosition!="0":
				self.tTempExceptedPosition=self.sExceptedPosition.split(self.sOptionInternalItemSeparator)
				try:
					self.tExceptedPosition=[int(X) for X in self.tTempExceptedPosition]
					#IsInt(tExceptedPosition,"position no take in count for the calcul of identity")
				except ValueError:
					exit("Error : ExceptedPosition must be defined by integer, process broken")
			else:
				self.tExceptedPosition=[]	
		
			self.iMatrixSizeWithoutExcepted=self.iMatrixSize-len(self.tExceptedPosition)
			self.fTheoricalIdentityMax=float(self.iMatrixSizeWithoutExcepted*2)
	
			if self.fTheoricalIdentityMax < self.fIdentityMax:
				print("Warning : HomologyMax can't be superior than 2x matrix length, default value "+str(self.fTheoricalIdentityMax)+" used")
				self.fIdentityMax=self.fTheoricalIdentityMax

	def result_file_name(self) :
		#generateur du nom de fichier resultat
		sDisplayData="NoDisplaySp"
		if len(self.tDisplaySpecies)!=0:
			sDisplayData="DisplaySp"+"-".join(self.tDisplaySpecies)
		sOrthoData="NoOrthoSp"
		if len(self.tOrthoSpecies)!=0:
			sOrthoData="OrthoSp"+"-".join(self.tOrthoSpecies)
		sMatchData="NoMatchSp"
		if len(self.tMatchSpecies)!=0:
			sMatchData="MatchSp"+"-".join(self.tMatchSpecies)
		sRefData="refSp"+self.sRefSpecies

		sValueData="ValueFilter"
		if self.fMinScore!=0:
			sValueData+="-Min"+self.sScoreValueMin
		if self.sScoreValueMax!="Nothing":
			sValueData+="-Max"+self.sScoreValueMax
		if sValueData=="ValueFilter":
			sValueData="NoValueFilter"

		sIdentityData="NoIdentityFilter"
		if self.bCheckHomogeneous:
			sIdentityData="IdentityFilter-"+str(self.fIdentityMin)+"-"+str(self.fIdentityMax)

		if self.bSingleHitAnalysis:
			sIdentityData+="-BestHit"
		else:
			sIdentityData+="-AllHit"
		if len(self.tExceptedPosition)!=0:
			sIdentityData+="-ExceptedPosition-"+"-".join(tTempExceptedPosition)
		else:
			sIdentityData+="-NoExceptedPosition"

		#sResultFileName="CompaMatrix.tsv"	
		self.sResultFileName="CompaMatrix-"+sVersionScript+"-ForJob"+self.sJobId+"_"+sValueData+"_"+self.sRegion+"_"+sRefData+"_"+sMatchData+"_"+sOrthoData+"_"+sDisplayData+"_"+sIdentityData

		print("default result file name used : "+self.sResultFileName)


class SQLquery :

	def __init__(self, cmd, arg) :
		self.result = SubmitMySQLCommand(cmd,arg.sHost,arg.sDatabase,arg.sUser,arg.sPw)


class Analysis :
	def __init__(self,arg) :
		self.dGenome2Name = {}
		self.dClique2HGNC = {}
		self.dClique2Biotype = {}
		self.dGene2Data = {}
		self.dGenome2GeneSample = {}
		self.dCheckInclude = {}
		self.dLost = {}
		self.dClique2Gene = {}
		self.dDisplay2Ortholog = {}
		self.dGenome2ENStag = {}	
		self.dClique2UniqPvalue = {}
		self.dFinal = {}
		#Take in count the tool_id (target the good column)
		sSQLcommand="select Tag from "+arg.sVersionDb+"_job j, "+arg.sVersionDb+"_tool t where j.tool_id=t.tool_id and j.job_id="+arg.sJobId+";"
		#tSQLtable=SubmitMySQLCommand(sSQLcommand,arg.sHost,arg.sDatabase,arg.sUser,arg.sPw)
		query = SQLquery(sSQLcommand,arg)
		#sToolColTarget=tSQLtable[0][0]
		self.sToolColTarget = query.result[0][0]

	def getAllSpecies(self,arg) : #step1
		sSQLcommand="select genome_db_id,name from "+arg.sVersionDb+"_genome_db;"
		#tSQLtable=SubmitMySQLCommand(sSQLcommand,sHost,sDatabase,sUser,sPw)
		query = SQLquery(sSQLcommand,arg)
		for tRows in query.result :
			sGenomeId=tRows[0]
			sSpeciesName=tRows[1]
			if sGenomeId in arg.tAllSpecies:
				self.dGenome2Name[sGenomeId]=sSpeciesName

	def getHGNCAnnotations(self,arg) : #step2
		sSQLcommand="select clique_id, fullName, shortName from "+arg.sVersionDb+"_HGNC h, "+arg.sVersionDb+"_asso_HGNC_clique a where h.HGNC_id=a.HGNC_id;"
		query = SQLquery(sSQLcommand,arg)
		for tRows in query.result :
			sCliqueId=tRows[0]
			sFullName=" ".join(tRows[1:-1])
			sShortName=tRows[-1]
			self.dClique2HGNC[sCliqueId]={"full":sFullName,"short":sShortName}

	def getCliqueBiotypeAnnotation(self,arg) : #step3
		sSQLcommand="select c.clique_id, t.biotype from "+arg.sVersionDb+"_gene_"+arg.sRefSpecies+" t, "+arg.sVersionDb+"_gene g, "+arg.sVersionDb+"_clique_element c where t.stable_id = g.stable_id and g.gene_id = c.sp"+arg.sRefSpecies+";"
		query = SQLquery(sSQLcommand,arg)
		for tRows in query.result :
			sCliqueId=tRows[0]
			sBiotype=tRows[1]
			self.dClique2Biotype[sCliqueId]=sBiotype

	def getRefHits(self,arg) :  #step4, merge with other species ?

		sSQLcommand="select gene_id,start_match,end_match,sequence_match,strand_match,"+self.sToolColTarget+",scan_id from "+arg.sVersionDb+"_scan_"+arg.sRefSpecies+" s where job_id="+arg.sJobId+" "+arg.sregionClause+" and s.match=1;"
		query = SQLquery(sSQLcommand,arg)
		for tRows in query.result :
			sGeneId=tRows[0]
			dData={"Start":tRows[1],"End":tRows[2],"Sequence":tRows[3],"Strand":tRows[4],self.sToolColTarget:tRows[5]}
			sScanId=tRows[6]
			fCurrentScore=float(tRows[5])
			bAddHit=True

			if arg.bApplyScoreFilter:
				if arg.bApplyMinScore:
					if fCurrentScore < arg.fMinScore:
						bAddHit=False
				if arg.bApplyMaxScore:
					if fCurrentScore > arg.fMaxScore:
						bAddHit=False
			if bAddHit:
				if arg.bSingleHitAnalysis:
					try:
						fStoredScore = float(self.dGene2Data[sGeneId][self.dGene2Data[sGeneId].keys()[0]][self.sToolColTarget])
						if fCurrentScore < fStoredScore:
							self.dGene2Data[sGeneId]={sScanId:dData}
					except KeyError:
						self.dGene2Data[sGeneId]={sScanId:dData}
				else:
					try:
						self.dGene2Data[sGeneId][sScanId]=dData
					except KeyError:
						self.dGene2Data[sGeneId]={sScanId:dData}

		sSpeciesId=sGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId             
 		self.dGenome2GeneSample[sSpeciesId]=sGeneId


	def getRefCliques(self,arg) : #step 5
		tAvailableColumn= ["sp"+sId for sId in arg.tAllSpecies]
		sSQLcommand="select clique_id,"+",".join(tAvailableColumn)+" from "+arg.sVersionDb+"_clique_element where sp"+arg.sRefSpecies+" is not null;"
		query = SQLquery(sSQLcommand, arg)
		for tRows in query.result :
			sCliqueId=tRows[0]
			sRefGeneId=tRows[1]
			tListOfOtherGene=[geneId for geneId in tRows[2:] if geneId!="NULL"]
			try:
				dTempVar=self.dGene2Data[sRefGeneId]
			except KeyError:
				#No valide refGene, pass
				continue
			bFirstApparition=False
			try:
				dTempVar=self.dCheckInclude[sRefGeneId]
			except KeyError:
				bFirstApparition=True
			if bFirstApparition:
				self.dCheckInclude[sRefGeneId]={sCliqueId:tListOfOtherGene}
			else:
				bHaveChange=True
				while bHaveChange:
					bIncludeOrSimilar=False
					bRemoveCurrentClique=False
					bRemoveStoredClique=False

					for sStoredClique in self.dCheckInclude[sRefGeneId]:
						bIncludeOrSimilar=False
						bRemoveCurrentClique=False
						bRemoveStoredClique=False

						tStoredListOfOtherGene=self.dCheckInclude[sRefGeneId][sStoredClique]
						tGeneSpecificToStoredClique=GetDifferenceBetweenList(tStoredListOfOtherGene,tListOfOtherGene)
						tGeneSpecificToCurrentClique=GetDifferenceBetweenList(tListOfOtherGene,tStoredListOfOtherGene)

						if len(tGeneSpecificToStoredClique)==0 and len(tGeneSpecificToCurrentClique)==0:
							bIncludeOrSimilar=True
							bRemoveCurrentClique=True
						elif len(tGeneSpecificToStoredClique)==0:
							bIncludeOrSimilar=True
							bRemoveStoredClique=True
						elif len(tGeneSpecificToCurrentClique)==0:
							bIncludeOrSimilar=True
							bRemoveCurrentClique=True
						else:
							bIncludeOrSimilar=False
						if bIncludeOrSimilar:
							break

					if not bRemoveCurrentClique and not bRemoveStoredClique:
						self.dCheckInclude[sRefGeneId][sCliqueId]=list(tListOfOtherGene)
						bHaveChange=False
					elif bRemoveStoredClique:
 						del self.dCheckInclude[sRefGeneId][sStoredClique]
 					elif bRemoveCurrentClique:
 						bHaveChange=False
				
	def orthoCliquesFilter(self,arg) : #step6
		for sRefGeneId in self.dCheckInclude:
			for sCliqueId in self.dCheckInclude[sRefGeneId]:
				tListOfGenes = self.dCheckInclude[sRefGeneId][sCliqueId]
				tListOfRepresentedSpecies=[]
				for sSpeciesGeneId in tListOfGenes:
					sSpeciesId = sSpeciesGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
					tListOfRepresentedSpecies.append(sSpeciesId)

				bValidClique=True
				for sMatchSpecies in arg.tMatchSpecies:
					if sMatchSpecies not in tListOfRepresentedSpecies:
						bValidClique=False
						try :
							self.dLost[sMatchSpecies]+=1
						except KeyError:
							self.dLost[sMatchSpecies]=1
				for sOrthoSpecies in arg.tOrthoSpecies:
					if sOrthoSpecies not in tListOfRepresentedSpecies:
						bValidClique=False
						try:
							self.dLost[sOrthoSpecies]+=1
						except KeyError:
							self.dLost[sOrthoSpecies]=1
				for sDisplaySpecies in arg.tDisplaySpecies:
					if not sDisplaySpecies in self.dDisplay2Ortholog :
						self.dDisplay2Ortholog[sDisplaySpecies]={}
					if sDisplaySpecies in tListOfRepresentedSpecies:
						self.dDisplay2Ortholog[sDisplaySpecies][sCliqueId]="one2one"
					else:
						self.dDisplay2Ortholog[sDisplaySpecies][sCliqueId]="no ortholog"
				if bValidClique:
					tNewListOfGenes=[sRefGeneId]+tListOfGenes
					tNewListOfGenes=sorted(tNewListOfGenes)
					dbListOfGenes=tuple(tNewListOfGenes)
					self.dClique2Gene[sCliqueId]=dbListOfGenes

	def getAllHits(self,arg) : #step 7, merge with step 4 ?	
		iCount=1
		sOneTarget = arg.sAbsentValue
		for sSpecies in self.dGenome2Name.keys():
			if sSpecies== arg.sRefSpecies:
				#Already processed
				continue
			print str(iCount)+"/"+str(len(self.dGenome2Name)-1)
			iCount+=1
			sSQLcommand="select gene_id,start_match,end_match,sequence_match,strand_match,"+self.sToolColTarget+",scan_id from "+arg.sVersionDb+"_scan_"+sSpecies+" s where job_id="+arg.sJobId+" "+arg.sregionClause+" and s.match=1;"
			query = SQLquery(sSQLcommand, arg)
			for tRows in query.result:
				sGeneId=tRows[0]
				sSpeciesId=sGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
				dData={"Start":tRows[1],"End":tRows[2],"Sequence":tRows[3],"Strand":tRows[4],self.sToolColTarget:tRows[5]}
				sScanId=tRows[6]
				fCurrentScore=float(tRows[5])

				bAddHit=True

				if arg.bApplyScoreFilter:
					#if sSpeciesId in tDisplayToFilter or not sSpeciesId in tDisplaySpecies: 
					if not sSpeciesId in arg.tDisplaySpecies:
						if arg.bApplyMinScore:
							if fCurrentScore < arg.fMinScore:
								bAddHit=False
						if arg.bApplyMaxScore:
							if fCurrentScore > arg.fMaxScore:
								bAddHit=False

				if bAddHit :
					if arg.bSingleHitAnalysis:
							try:
								fStoredScore=float(self.dGene2Data[sGeneId][self.dGene2Data[sGeneId].keys()[0]][self.sToolColTarget])
								if fCurrentScore < fStoredScore:
									self.dGene2Data[sGeneId]={sScanId:dData}
							except KeyError:
								self.dGene2Data[sGeneId]={sScanId:dData}
					else:
						try:
							self.dGene2Data[sGeneId][sScanId]=dData
						except KeyError:
							self.dGene2Data[sGeneId]={sScanId:dData}

			sSpeciesId=sGeneId[:-11]              
			self.dGenome2GeneSample[sSpeciesId]=sGeneId
			self.sOneTarget=tRows[3]

	def hitCliquesFilter(self,arg) : #step 8
		tUnalyzedGenes=[]
		if self.sOneTarget == arg.sAbsentValue:
 			print("Warning : no valide results, step skipped")
 		elif not arg.bGetResult:
			print("Warning : Count option activated, no homologoy tested, step skipped")
			self.dClique2Score={}
		elif bCheckHomogeneous:
			self.dClique2Score={}
			iCount=1
			iFinalCount=len(dClique2Gene)
			for sCliqueId in dClique2Gene.keys():
				if iCount%10==0:
					print str(iCount)+"/"+str(iFinalCount)
				iCount+=1
				print sCliqueId
 				try:
					tTargetedGenes=[X for X in dClique2Gene[sCliqueId] if X[:-11] in tCheckedSpecies]
					self.dWorkGene2Data={}
					for sGeneId in tTargetedGenes:
						self.dWorkGene2Data[sGeneId]=dGene2Data[sGeneId]

				except KeyError:
					#one or more species have no valid hit for the clique. Clique is removed
					del dClique2Gene[sCliqueId]
					try:
						dLost["Heterogeneous"]+=1
					except KeyError:
						dLost["Heterogeneous"]=1
					continue
				#Remove Clique if too much combination available
				iCombination=1
				bBreak=False
				for sGeneId in self.dWorkGene2Data:
					#print(dWorkGene2Data[sGeneId])
					iCombination=iCombination*len(dWorkGene2Data[sGeneId])
					if iCombination>iHomogeneityLimit:
						print(sGeneId+" is rejected : too much hits to compute homogenity ("+str(iCombination)+")")
						bBreak=True
						break
				if bBreak:
					try:
						tUnalyzedGenes.append(dClique2HGNC[sCliqueId]["short"])
					except KeyError:
						tUnalyzedGenes.append("NoHGNC OrthocisId-"+sGeneId)
					#print len(tUnalyzedGenes)
					del dClique2Gene[sCliqueId]
					try:
						dLost["Heterogeneous"]+=1
					except KeyError:
						dLost["Heterogeneous"]=1
					continue
        #

					dbTemp=launchHomology4MultiHits(HOMOLOGY_SCRIPT_MODE,self.dWorkGene2Data,tExceptedPosition,fIdentityMax)

					fIdentityScore=dbTemp[0]
					dBestHit4Genes=dbTemp[1]

					if fIdentityScore>fIdentityMax or fIdentityScore<fIdentityMin:
						del dClique2Gene[sCliqueId]
						try:
							dLost["Heterogeneous"]+=1
						except KeyError:
							dLost["Heterogeneous"]=1
						continue

					dClique2Score[sCliqueId]=fIdentityScore
					for sGeneId in dBestHit4Genes:
						dGene2Data[sGeneId]=dBestHit4Genes[sGeneId]
		else:
			print("No identityThreshold fixed, step skipped")

	def getEnsemblId(self,arg) : #step9
	        argid = self.dGenome2GeneSample.values()
        	if argid is not None and len(argid)>0:
			sSQLcommand="select stable_id,genome_db_id from "+arg.sVersionDb+"_gene where gene_id="+" or gene_id=".join(self.dGenome2GeneSample.values())+";"
			query = SQLquery(sSQLcommand,arg)
			for tRows in query.result :
				sENSid=tRows[0]
				sENStag=sENSid[:-11]
				sGenomeId=tRows[1]
				self.dGenome2ENStag[sGenomeId]=sENStag

	def getClique2Pvalue(self,arg) : #step10 
		for sCliqueId in self.dClique2Gene:
			fWorstGlobalPvalue=-1
			for sGeneId in self.dClique2Gene[sCliqueId]:
				if sGeneId[:-11] in arg.tCheckedSpecies:
					#Warning: a gene can have no hit, so no occurences in dGene2Data
					try:
 						fLocalBestPvalue=1.0
						for sScanId in self.dGene2Data[sGeneId]:
							fGenePvalue=float(self.dGene2Data[sGeneId][sScanId][self.sToolColTarget])
							fLocalBestPvalue=min(fGenePvalue,fLocalBestPvalue)
						fWorstGlobalPvalue=max(fLocalBestPvalue,fWorstGlobalPvalue)
					except KeyError:
						continue
			self.dClique2UniqPvalue[sCliqueId]=fWorstGlobalPvalue


	def getFinalDict(self,arg) : #step11
		iCount=1
		for sCliqueId in self.dClique2Gene:
			if iCount%100==0:
				print str(iCount)+"/"+str(len(self.dClique2Gene))
			iCount+=1
			#Warning: a clique can have no HGNC annotation
			try:
				self.dFinal[sCliqueId]={"HGNC":{"full":self.dClique2HGNC[sCliqueId]["full"],"short":self.dClique2HGNC[sCliqueId]["short"]},"Content":{}}
			except KeyError:
				self.dFinal[sCliqueId]={"HGNC":{"full":"no HGNC","short":"no HGNC"},"Content":{}}
			if arg.bCheckHomogeneous:
				self.dFinal[sCliqueId]["Score"] = self.dClique2Score[sCliqueId]
			self.dFinal[sCliqueId]["Pvalue"] = self.dClique2UniqPvalue[sCliqueId]
			bValidClique=True
			bRefFailure=False
			tMatchFailure=[]
			for sGeneId in self.dClique2Gene[sCliqueId]:
				#Warning: a gene can have no hit, so no occurences in dGene2Data
				self.dFinal[sCliqueId]["Content"][sGeneId]={}
				try:
					for sScanId in self.dGene2Data[sGeneId]:
						self.dFinal[sCliqueId]["Content"][sGeneId][sScanId] = self.dGene2Data[sGeneId][sScanId]
				except KeyError:
					sSpeciesId=sGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
					if sSpeciesId == arg.sRefSpecies:
						bValidClique=False
						bRefFailure=True
					elif sSpeciesId in arg.tMatchSpecies:
						bValidClique=False
						tMatchFailure.append(sSpeciesId)
					else:
						self.dFinal[sCliqueId]["Content"][sGeneId]["noHit"]={"Start":"","End":"","Sequence":"","Strand":"",self.sToolColTarget:"NA"}
			if not bValidClique:
				del self.dFinal[sCliqueId]
				if not bRefFailure:
					for sSpeciesId in tMatchFailure:
						try:
							self.dLost[sSpeciesId]+=1
						except KeyError:
							self.dLost[sSpeciesId]=1

	def getConservationValue(self,arg) : #step12
		if arg.bSingleHitAnalysis or arg.bCheckHomogeneous:
			print "step 12"
			for sCliqueId in self.dFinal:
				self.dFinal[sCliqueId]["ConservationValue"]=GetCommonValue(self.dFinal[sCliqueId]["Content"],tExceptedPosition,tCheckedSpecies)
 		else:
			print("Warning : many hits conserved by gene, no ConservationValue computed")

	def countWrite(self,arg,finalRes) : #step 13 count mode
		#result = count.json
		if self.dLost is None:
			self.dLost=dict()
		if "Heterogeneous" not in self.dLost.keys():
			self.dLost["Heterogeneous"]=None
		finalRes.dSpeciesCount={"Global":{"cliqueNumber":0,"geneNumber":0,"hitNumber":0,"heterogeneousRejected":self.dLost["Heterogeneous"]},"Reference":{"geneNumber":0,"hitNumber":0},"MatchSpecies":{},"OrthoSpecies":{},"DisplaySpecies":{}}
		dSpeciesClassification={}
		dSpeciesClassification[arg.sRefSpecies]={"type":"Reference","name":self.dGenome2Name[arg.sRefSpecies]}
		for sSpeciesId in arg.tMatchSpecies:
			dSpeciesClassification[sSpeciesId]={"type":"MatchSpecies","name":self.dGenome2Name[sSpeciesId]}
			try:
				iLostNumber=self.dLost[sSpeciesId]
			except KeyError:
				iLostNumber=0
			finalRes.dSpeciesCount["MatchSpecies"][self.dGenome2Name[sSpeciesId]]={"geneNumber":0,"hitNumber":0,"lostClique":iLostNumber}
		for sSpeciesId in arg.tOrthoSpecies:
			dSpeciesClassification[sSpeciesId]={"type":"OrthoSpecies","name":self.dGenome2Name[arg.sRefSpecies]}
			try:
				iLostNumber=self.dLost[sSpeciesId]
			except KeyError:
				iLostNumber=0
			finalRes.dSpeciesCount["OrthoSpecies"][self.dGenome2Name[sSpeciesId]]={"geneNumber":0,"hitNumber":0,"lostClique":iLostNumber}
		for sSpeciesId in arg.tDisplaySpecies:
			dSpeciesClassification[sSpeciesId]={"type":"DisplaySpecies","name":self.dGenome2Name[arg.sRefSpecies]}
			finalRes.dSpeciesCount["DisplaySpecies"][self.dGenome2Name[sSpeciesId]]={"geneNumber":0,"hitNumber":0}

		for sCliqueId in self.dFinal:
			finalRes.dSpeciesCount["Global"]["cliqueNumber"]+=1
			for sGeneId in self.dFinal[sCliqueId]["Content"]:
				finalRes.dSpeciesCount["Global"]["geneNumber"]+=1
				sSpeciesId=sGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
				sSpeciesName=self.dGenome2Name[sSpeciesId]
				iHitNumber=len(self.dFinal[sCliqueId]["Content"][sGeneId])
				finalRes.dSpeciesCount["Global"]["hitNumber"]+=iHitNumber
				if sSpeciesId == arg.sRefSpecies:
					finalRes.dSpeciesCount["Reference"]["geneNumber"]+=1
					finalRes.dSpeciesCount["Reference"]["hitNumber"]+=iHitNumber
					finalRes.internal_genes.append(sGeneId)
				if sSpeciesId in arg.tMatchSpecies:
					finalRes.dSpeciesCount["MatchSpecies"][sSpeciesName]["geneNumber"]+=1
					finalRes.dSpeciesCount["MatchSpecies"][sSpeciesName]["hitNumber"]+=iHitNumber
				if sSpeciesId in arg.tOrthoSpecies:
					finalRes.dSpeciesCount["OrthoSpecies"][sSpeciesName]["geneNumber"]+=1
					finalRes.dSpeciesCount["OrthoSpecies"][sSpeciesName]["hitNumber"]+=iHitNumber
				if sSpeciesId in arg.tDisplaySpecies:
					finalRes.dSpeciesCount["DisplaySpecies"][sSpeciesName]["geneNumber"]+=1
					finalRes.dSpeciesCount["DisplaySpecies"][sSpeciesName]["hitNumber"]+=iHitNumber
		finalRes.computeGenes(200, arg)
		finalRes.dSpeciesCount["Reference"]["selected_genes"]=finalRes.genes
		print("sJsonOutputName"+arg.sJsonOutputName)
		FILE=open(arg.sJsonOutputName,"w")
		FILE.write(json.dumps(finalRes.dSpeciesCount, indent=4))
		FILE.close()

	def resultWrite(self,arg) :#step13 result mode	
		#Add BestPvalue column
		bBestPvalue=False
		if not arg.bSingleHitAnalysis and not arg.bCheckHomogeneous:
			bBestPvalue=True
		#Create headline of the file
		sHeaderLine="HGNC long\tHGNC short\tBiotype\t"
		if len(tExceptedPosition)==0:
			sExceptedTag=""
		else:
			sExceptedTag="[excepted-"+"-".join(arg.tExceptedPosition)+"]"
		if arg.bCheckHomogeneous:
			sHeaderLine+="HomologyScore"+sExceptedTag+"\t"
		if arg.bCheckHomogeneous or arg.bSingleHitAnalysis:
			sHeaderLine+="ConservationValue"+sExceptedTag+"\t"
		sHeaderLine+="WorstPvalue\t"
		for sSpeciesId in arg.tAllSpecies:
			sSpeciesShortId=""
			if sSpeciesId in self.dGenome2ENStag:
				sSpeciesShortId = self.dGenome2ENStag[sSpeciesId][:-1]
				sSpeciesShortId = sSpeciesShortId.replace("ENS","")
			if sSpeciesShortId=="":
				sSpeciesShortId="HSA"
			sSpeciesCategory = sSpeciesShortId+"."
			if sSpeciesId == arg.sRefSpecies:
				sSpeciesCategory+="Ref"
			elif sSpeciesId in arg.tMatchSpecies:
				sSpeciesCategory+="Match"
			elif sSpeciesId in arg.tOrthoSpecies:
				sSpeciesCategory+="Ortho"
			elif sSpeciesId in arg.tDisplaySpecies:
				sSpeciesCategory+="Display"
			else:
				exit("Unknow speciesId "+sSpeciesId+", process broken")
			sHeaderLine+=sSpeciesCategory
			if bBestPvalue:
				sHeaderLine+="\t"+sSpeciesShortId+".BestPvalue"
			sHeaderLine+="\t"+sSpeciesShortId+".gene\t"+sSpeciesShortId+"."+sToolColTarget+"\t"+sSpeciesShortId+".coord\t"+sSpeciesShortId+".sequence\t"
			if sSpeciesId in arg.tDisplaySpecies:
				sHeaderLine+=sSpeciesShortId+".ortholog\t"

		FILE=open(sResultFileName+".tsv","w")
		FILE.write(sHeaderLine+"\n")

		#Write each line, 1 by valide clique
		iCount=1
		iWritedClique=0
		for sCliqueId in self.dFinal:
			if iCount%100==0:
				print str(iCount)+"/"+str(len(self.dFinal))
			iCount+=1
			#HGNC annotation
			sCurrentLine = self.dFinal[sCliqueId]["HGNC"]["full"]+"\t"+self.dFinal[sCliqueId]["HGNC"]["short"]+"\t"
			#Biotype annotation
			sCurrentLine += self.dClique2Biotype[sCliqueId]+"\t"
			#Identity Score
			if bCheckHomogeneous:
				sCurrentLine+=str(round(self.dFinal[sCliqueId]["Score"],2))+"\t"
			if arg.bCheckHomogeneous or arg.bSingleHitAnalysis:
				sCurrentLine+=str(self.dFinal[sCliqueId]["ConservationValue"])+"\t"
			#Pvalue Score
			sCurrentLine += str(self.dFinal[sCliqueId]["Pvalue"])+"\t"
			#For each species
			for sSpeciesId in arg.tAllSpecies:
				asens=0
				if sSpeciesId in self.dGenome2ENStag:
					asens=1
				else:
					asens=0
					self.dGenome2ENStag[sSpeciesId]=""
				#CategorySpacer
				sCurrentLine+="\t"
				bNoGeneValide=True
				#Search the affiliate gene
				for sGeneId in dFinal[sCliqueId]["Content"]:
					if sGeneId[:-11]!=sSpeciesId: #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
						continue
					bNoGeneValide=False
					#two cases : have at least one valide hit or have no hit
					#but if it is a display or an Ortho, we must indicate the gene name
					if dFinal[sCliqueId]["Content"][sGeneId].keys()==["noHit"]:
						#if displaySpecies, add ortholog information
						sCurrentLine+=dGenome2ENStag[sSpeciesId]+sGeneId[-11:]+"\t\t\t\t"
						if bBestPvalue:
							sCurrentLine+="\t"
						if sSpeciesId in tDisplaySpecies:
							sCurrentLine+=dDisplay2Ortholog[sSpeciesId][sCliqueId]+"\t"
					else:
						sENSid=dGenome2ENStag[sSpeciesId]+sGeneId[-11:]
						fBestPvalue=1.0
						sScore=""
						sCoord=""
						sSeq=""
						dElement={}
						for sScanId in dFinal[sCliqueId]["Content"][sGeneId]:
							fBestPvalue=min(float(dGene2Data[sGeneId][sScanId][sToolColTarget]),fBestPvalue)
							try:
								dElement[float(dFinal[sCliqueId]["Content"][sGeneId][sScanId][sToolColTarget])]=(dFinal[sCliqueId]["Content"][sGeneId][sScanId]["Start"]+":"+dFinal[sCliqueId]["Content"][sGeneId][sScanId]["End"]+"("+dFinal[sCliqueId]["Content"][sGeneId][sScanId]["Strand"]+")",dFinal[sCliqueId]["Content"][sGeneId][sScanId]["Sequence"])
							except ValueError:
								print dFinal[sCliqueId]["Content"][sGeneId]
								exit()
						bFirst=True
						for fScore in sorted(dElement.keys()):
							if bFirst:
								bFirst=False
								sSep=""
							else:
								sSep=" "
							sScore+=sSep+str(fScore)
							sCoord+=sSep+dElement[fScore][0]
							sSeq+=sSep+dElement[fScore][1]
						sCurrentLine+=sENSid
						if bBestPvalue:
							sCurrentLine+="\t"+str(fBestPvalue)
						sCurrentLine+="\t"+sScore+"\t"+sCoord+"\t"+sSeq+"\t"
						if sSpeciesId in tDisplaySpecies:
							sCurrentLine+=dDisplay2Ortholog[sSpeciesId][sCliqueId]+"\t"
                                #end of the loop
				if bNoGeneValide:
					sCurrentLine+="\t\t\t\t"
					if sSpeciesId in tDisplaySpecies:
						sCurrentLine+=dDisplay2Ortholog[sSpeciesId][sCliqueId]+"\t"
			#End of the line                
			sCurrentLine+="\n"
			FILE.write(sCurrentLine)
			iWritedClique+=1

		FILE.close()

		print sResultFileName
		print("The file "+sResultFileName+" contains "+str(iWritedClique)+" genes")
		if len(tUnalyzedGenes)>0:
			sUnalyzedGenes="UnalyzedGenes-lim"+str(iHomogeneityLimit)+"-"+sResultFileName+".txt"
			FILE=open(sUnalyzedGenes,"w")
			FILE.write("\n".join(tUnalyzedGenes))
			FILE.close()
			print("There is "+str(len(tUnalyzedGenes))+" genes with too much hits for the homogenity computation (limit: "+str(iHomogeneityLimit)+").")
			print("This genes are stored in the file "+sUnalyzedGenes)

				
		
def mainprocessclasses() :

	#Step0:Load all arguments

	arg = Arguments() #the arg object will contain all arguments
	arg.add_constants() #add constants (they never change)
	arg.which_varsource(arg.options) #which var source is used ? json file, json path or json string ?
	arg.add_arguments(arg.source) #arguments loaded from previous json data
	arg.add_database_parameters() #load SQL parameters
	arg.check_validity() #check if arguments are valid
	arg.concatenate() #concatenate all arguments
	arg.result_file_name() #create output_file_name (for result mode)

	print "Arguments loaded"

	finalRes = FinalDataStruct(arg)

	analysis = Analysis(arg)
	
	print "Analysis object created"

	iStepCount = 1

	#Step0:get matrix size sJsonKeyMatrixSize
	iStepCount=StepIncrease(iStepCount,"Get matrix size",0)
	arg.getMatrixSize()

	#Step1:retrieve all available species
	iStepCount = StepIncrease(iStepCount,"Retrieve all species",0)
	analysis.getAllSpecies(arg)

	#Step2:retrieve all HGNC annotation for clique
	iStepCount = StepIncrease(iStepCount,"Retrieve all HGNC",0)
	analysis.getHGNCAnnotations(arg)

	#Step3:retrieve all biotype annotation for clique
	iStepCount = StepIncrease(iStepCount,"Retrieve all biotype annotation for clique",0)
	analysis.getCliqueBiotypeAnnotation(arg)

	#Step4: retrieve all hits on gene for the jobId
	iStepCount = StepIncrease(iStepCount,"Retrieve all hits on ref species",0)
	analysis.getRefHits(arg)

	#Step5: retrieve all clique where sRefSpecies have a gene and remove clique that are included in other clique ac     cording the species selection
        iStepCount=StepIncrease(iStepCount,"Retrieve all cliques",0)
	analysis.getRefCliques(arg)

        #Step6: Check that ortho/matchSpecies have a gene in each clique
        #               and note if display have or not a ortholog gene
	iStepCount=StepIncrease(iStepCount,"Take in count Match/Ortho/Display specificity",0)
	analysis.orthoCliquesFilter(arg)

        #Step7:retrieve all hits on gene for the jobId
        #NB: if bSingleHitAnalysis is True, retrieve only the best hit
        iStepCount=StepIncrease(iStepCount,"Retrieve all hits for each species",0)
	analysis.getAllHits(arg)
	
        #Step8: Compute and check if hits are homogeneous for MatchSpecies and RefSpecies
        iStepCount=StepIncrease(iStepCount,"Apply homology filter",0)
	analysis.hitCliquesFilter(arg)

        #Step9:retriever the ENS id for each species
        iStepCount=StepIncrease(iStepCount,"Retrieve EnsemblId",0)
	analysis.getEnsemblId(arg)

        #Step10:set the global pvalue for each clique
        iStepCount=StepIncrease(iStepCount,"Set the global pvalue for each clique",0)
	analysis.getClique2Pvalue(arg)

        #Step11:create the final dictionnary
        iStepCount=StepIncrease(iStepCount,"Regroup all results",0)
	analysis.getFinalDict(arg)

        #Step12:
        iStepCount=StepIncrease(iStepCount,"Compute the ConservationValue",0)
	analysis.getConservationValue(arg)

        #Step13:write the result
        iStepCount=StepIncrease(iStepCount,"Write the results",0)
        if not arg.bGetResult :
		analysis.countWrite(arg,finalRes)
	else :
		analysis.resultWrite(arg)


mainprocessclasses()

########fix :cleaning temp data
cmd=["bash","-c", "rm -f ./TempCompa*" ]
ExecuteBashCommand(cmd)
cmd=["bash","-c", "rm -f ./datatmp*" ]
ExecuteBashCommand(cmd)
