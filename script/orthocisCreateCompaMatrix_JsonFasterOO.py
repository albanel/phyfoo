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
           	print sSQLcommand
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


class Clique : 

	def __init__(self, clique_id, refGene, otherGenes) :
		self.clique_id = clique_id
		self.refGene = refGene
		self.otherGenes = otherGenes

	#def add_refGenes_and_otherGenes(self, refGene, otherGenes) :
	#	self.refGene = refGene
	#	self.otherGenes = otherGenes

class Hit :

	def __init__(self, scan_id, gene_id, start, end, seq, strand, pvalue) :
		self.scan_id = scan_id
		self.gene_id = gene_id
		self.start = start
		self.end = end
		self.seq = seq
		self.strand = strand
		self.pvalue = float(pvalue)


class SQLquery :

	def __init__(self, cmd, oArgument) :
		self.result = SubmitMySQLCommand(cmd,oArgument.sHost,oArgument.sDatabase,oArgument.sUser,oArgument.sPw)


def mainprocess():

	#Load all arguments

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

	########################################################################
	#Take in count the tool_id (target the good column)
	sSQLcommand="select Tag from "+arg.sVersionDb+"_job j, "+arg.sVersionDb+"_tool t where j.tool_id=t.tool_id and j.job_id="+arg.sJobId+";"
	#tSQLtable=SubmitMySQLCommand(sSQLcommand,arg.sHost,arg.sDatabase,arg.sUser,arg.sPw)
	query = SQLquery(sSQLcommand,arg)
	#sToolColTarget=tSQLtable[0][0]
	sToolColTarget = query.result[0][0]
	
	'''
	########################################################################
	#DEBUG
	print "-----------DEBUG-----------"
	print "sSingleHitAnalysis:",sSingleHitAnalysis
	print "tExceptedPosition:",tExceptedPosition
	print "Homogeneous:",bCheckHomogeneous
	if bCheckHomogeneous: print "iMatrixSize:",iMatrixSize
	if bCheckHomogeneous: print "fIdentityMin:",fIdentityMin
	if bCheckHomogeneous: print "fIdentityMax:",fIdentityMax
	print "bApplyScoreFilter:",bApplyScoreFilter
	print "bApplyMinScore:",bApplyMinScore
	print "bApplyMaxScore:",bApplyMaxScore
	if bApplyMinScore: print "fMinScore:",fMinScore
	if bApplyMaxScore: print "fMaxScore:",fMaxScore
	print "-----------DEBUG-----------"
	'''
	########################################################################
	#Main
	iStepCount = 1
	#Step1:retrieve all available species
	iStepCount=StepIncrease(iStepCount,"Retrieve all species",0)
	dGenome2Name={}
	sSQLcommand="select genome_db_id,name from "+arg.sVersionDb+"_genome_db;"
	#tSQLtable=SubmitMySQLCommand(sSQLcommand,sHost,sDatabase,sUser,sPw)
	query = SQLquery(sSQLcommand,arg)
	for tRows in query.result :
		sGenomeId=tRows[0]
		sSpeciesName=tRows[1]

		if sGenomeId in arg.tAllSpecies:
			dGenome2Name[sGenomeId]=sSpeciesName

	#Step2:retrieve all HGNC annotation for clique
	iStepCount=StepIncrease(iStepCount,"Retrieve all HGNC",0)
	dClique2HGNC={}

        # warning : only for one2one (1 HNGC by clique)
	sSQLcommand="select clique_id, fullName, shortName from "+arg.sVersionDb+"_HGNC h, "+arg.sVersionDb+"_asso_HGNC_clique a where h.HGNC_id=a.HGNC_id;"
	#tSQLtable=SubmitMySQLCommand(sSQLcommand,sHost,sDatabase,sUser,sPw)
	query = SQLquery(sSQLcommand,arg)
	for tRows in query.result :
		sCliqueId=tRows[0]
		sFullName=" ".join(tRows[1:-1])
		sShortName=tRows[-1]

		dClique2HGNC[sCliqueId]={"full":sFullName,"short":sShortName}

	#Step3:retrieve all biotype annotation for clique
	iStepCount=StepIncrease(iStepCount,"Retrieve all biotype annotation for clique",0)
	dClique2Biotype={}
	sSQLcommand="select c.clique_id, t.biotype from "+arg.sVersionDb+"_gene_"+arg.sRefSpecies+" t, "+arg.sVersionDb+"_gene g, "+arg.sVersionDb+"_clique_element c where t.stable_id = g.stable_id and g.gene_id = c.sp"+arg.sRefSpecies+";"
	#tSQLtable=SubmitMySQLCommand(sSQLcommand,sHost,sDatabase,sUser,sPw)
	query = SQLquery(sSQLcommand,arg)
	for tRows in query.result :
		sCliqueId=tRows[0]
		sBiotype=tRows[1]

		dClique2Biotype[sCliqueId]=sBiotype

	
	#Step4: retrieve all hits on gene for the jobId
	#NB: if bSingleHitAnalysis is True, retrieve only the best hit
	iStepCount=StepIncrease(iStepCount,"Retrieve all hits on ref species",0)
	dGene2Data={}
	dGenome2GeneSample={}
	
	sSQLcommand="select gene_id,start_match,end_match,sequence_match,strand_match,"+sToolColTarget+",scan_id from "+arg.sVersionDb+"_scan_"+arg.sRefSpecies+" s where job_id="+arg.sJobId+" "+arg.sregionClause+" and s.match=1;"
	query = SQLquery(sSQLcommand,arg)
	for tRows in query.result :
		sGeneId=tRows[0]
		dData={"Start":tRows[1],"End":tRows[2],"Sequence":tRows[3],"Strand":tRows[4],sToolColTarget:tRows[5]}
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
		#bSingleHitAnalysis = True
		if bAddHit:
			if arg.bSingleHitAnalysis:
				try:
					fStoredScore = float(dGene2Data[sGeneId][dGene2Data[sGeneId].keys()[0]][sToolColTarget])
					if fCurrentScore < fStoredScore:
						dGene2Data[sGeneId]={sScanId:dData}
				except KeyError:
					dGene2Data[sGeneId]={sScanId:dData}
			else:
				try:
					dGene2Data[sGeneId][sScanId]=dData
				except KeyError:
					dGene2Data[sGeneId]={sScanId:dData}
	
		sSpeciesId=sGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId		
		dGenome2GeneSample[sSpeciesId]=sGeneId #useful only for result mode

	#print len(set(dGene2Data.keys()))
	'''
	#Step4 : OO version

	iStepCount=StepIncrease(iStepCount,"Retrieve all hits on ref species",0)
#	dGene2Data={}
	dGenome2GeneSample={}
	#sSQLcommand="select gene_id,start_match,end_match,sequence_match,strand_match,"+sToolColTarget+",scan_id from "+sVersionDb+"_scan_"+sRefSpecies+" s where job_id="+sJobId+" "+sregionClause+" and s.match=1;"
	
	#tSQLtable=SubmitMySQLCommand(sSQLcommand,sHost,sDatabase,sUser,sPw)

	sSQLcommand = "select gene_id,start_match,end_match,sequence_match,strand_match,"+sToolColTarget+",scan_id from "+arg.sVersionDb+"_scan_"+arg.sRefSpecies+" s where job_id="+arg.sJobId+" "+arg.sregionClause+" and s.match=1;"	

	query = SQLquery(sSQLcommand,arg)

	HitDict = {}

	for tRows in query.result :
		pvalue=float(tRows[5])
		sGeneId=tRows[0]		
		bAddHit=True

		if arg.bApplyScoreFilter :
			if arg.bApplyMinScore :
				if pvalue < arg.fMinScore :
					bAddHit=False
			if arg.bApplyMaxScore :
				if pvalue > arg.fMaxScore :
					bAddHit=False
		#print bAddHit
		#bSingleHitAnalysis=True
		if bAddHit :
			sScan_id=tRows[6]
			start=tRows[1]
			end=tRows[2]
			seq=tRows[3]
			strand=tRows[4]
			hitName=sScan_id+"_"+sGeneId+"_"+start+"_"+end+"_"+seq+"_"+strand+"_"+str(pvalue)
			hitName=Hit(sScan_id,sGeneId,start,end,seq,strand,pvalue) 
			if arg.bSingleHitAnalysis :
				try :
					fStoredScore=HitDict[sGeneId][0].pvalue #score of the lone stored hit
					if hitName.pvalue < fStoredScore :
						HitDict[sGeneId] = [hitName]
				except KeyError :
					HitDict[sGeneId] = [hitName]					
	
			else :
				#print "all hits kept"
				try : 
					HitDict[sGeneId].append(hitName)								
				except KeyError :
					HitDict[sGeneId]=[hitName]

		sSpeciesId=sGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId  
		dGenome2GeneSample[sSpeciesId]=sGeneId #useful only for result mode
	
	#print len(HitDict.keys())
	'''		
	#Step5: retrieve all clique where sRefSpecies have a gene and remove clique that are included in other clique according the species selection
	#Note:"Try-except" on dictionnary is faster than "if X in [data]-else" in case where len[data]>10000...
	iStepCount=StepIncrease(iStepCount,"Retrieve all cliques",0)
	tAvailableColumn= ["sp"+sId for sId in arg.tAllSpecies]
	dCheckInclude={}
	sSQLcommand="select clique_id,"+",".join(tAvailableColumn)+" from "+arg.sVersionDb+"_clique_element where sp"+arg.sRefSpecies+" is not null;"
	query = SQLquery(sSQLcommand, arg)	
	
	for tRows in query.result :
		sCliqueId=tRows[0]
		sRefGeneId=tRows[1]
		tListOfOtherGene=[geneId for geneId in tRows[2:] if geneId!="NULL"]
		try:
			dTempVar=dGene2Data[sRefGeneId]
		except KeyError:
			#No valide refGene, pass
			continue
		bFirstApparition=False
		try:
			dTempVar=dCheckInclude[sRefGeneId]
		except KeyError:
			bFirstApparition=True 
		if bFirstApparition:	
			dCheckInclude[sRefGeneId]={sCliqueId:tListOfOtherGene}
		else:
			bHaveChange=True
			while bHaveChange:
				bIncludeOrSimilar=False
				bRemoveCurrentClique=False
				bRemoveStoredClique=False
			
				for sStoredClique in dCheckInclude[sRefGeneId]:
					bIncludeOrSimilar=False
					bRemoveCurrentClique=False
					bRemoveStoredClique=False
				
					tStoredListOfOtherGene=dCheckInclude[sRefGeneId][sStoredClique]
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
					dCheckInclude[sRefGeneId][sCliqueId]=list(tListOfOtherGene)
					bHaveChange=False
				elif bRemoveStoredClique:
					del dCheckInclude[sRefGeneId][sStoredClique]
				elif bRemoveCurrentClique:
					bHaveChange=False

	'''	
	#Step 5 : OO version
	#ajout de tlist of other genes a la classe clique
	# try : dTempVar = dHitDict[gene_id]
	# si firstapparition : sClique_id=Clique(sClique_id,refgene,listofothergenes), clique_list.append(sClique_id)
	# sinon : for sStoredClique in clique_list : tStoredListOfOtherGene=sStoredClique.ListOfOtherGenes, reste est le meme
	iStepCount=StepIncrease(iStepCount,"Retrieve all cliques",0)
	tAvailableColumn= ["sp"+sId for sId in arg.tAllSpecies]
	dCheckInclude={}
	sSQLcommand="select clique_id,"+",".join(tAvailableColumn)+" from "+arg.sVersionDb+"_clique_element where sp"+arg.sRefSpecies+" is not null;"
	#tSQLtable=SubmitMySQLCommand(sSQLcommand,sHost,sDatabase,sUser,sPw)
	query = SQLquery(sSQLcommand,arg)	
	for tRows in query.result :
		sCliqueId=tRows[0]
		sRefGeneId=tRows[1]
		tListOfOtherGene=[geneId for geneId in tRows[2:] if geneId!="NULL"]
		try:
			dTempVar=HitDict[sRefGeneId]
		except KeyError:
			#No valide refGene, pass
			continue
		bFirstApparition=False
		try:
			dTempVar=dCheckInclude[sRefGeneId]
		except KeyError:
			bFirstApparition=True 
		sCliqueId=Clique(sCliqueId, sRefGeneId, tListOfOtherGene)
		if bFirstApparition:	
			#dCheckInclude[sRefGeneId]={sCliqueId:tListOfOtherGene}
			dCheckInclude[sRefGeneId]=[sCliqueId]
		else:
			bHaveChange=True
			while bHaveChange:
				bIncludeOrSimilar=False
				bRemoveCurrentClique=False
				bRemoveStoredClique=False
			
				for sStoredClique in dCheckInclude[sRefGeneId] :
					bIncludeOrSimilar=False
					bRemoveCurrentClique=False
					bRemoveStoredClique=False
				
					#tStoredListOfOtherGene=dCheckInclude[sRefGeneId][sStoredClique]
					tStoredListOfOtherGene=sStoredClique.otherGenes
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
					#dCheckInclude[sRefGeneId][sCliqueId]=list(tListOfOtherGene)
					dCheckInclude[sRefGeneId].append(sCliqueId)
					bHaveChange=False
				elif bRemoveStoredClique:
					#del dCheckInclude[sRefGeneId][sStoredClique]
					dCheckInclude[sRefGeneId].remove(sStoredClique)
				elif bRemoveCurrentClique:
					bHaveChange=False

	'''		
	#Step6: Check that ortho/matchSpecies have a gene in each clique
	#		and note if display have or not a ortholog gene
	iStepCount=StepIncrease(iStepCount,"Take in count Match/Ortho/Display specificity",0)
	dLost={}
	dClique2Gene={}
	dDisplay2Ortholog={}
	for sRefGeneId in dCheckInclude:
		for sCliqueId in dCheckInclude[sRefGeneId]:
			tListOfGenes=dCheckInclude[sRefGeneId][sCliqueId]
			tListOfRepresentedSpecies=[]
			for sSpeciesGeneId in tListOfGenes:
				sSpeciesId=sSpeciesGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
				tListOfRepresentedSpecies.append(sSpeciesId)
		
			bValidClique=True
			for sMatchSpecies in arg.tMatchSpecies:
				if sMatchSpecies not in tListOfRepresentedSpecies:
					bValidClique=False
					try:
						dLost[sMatchSpecies]+=1
					except KeyError:
						dLost[sMatchSpecies]=1
			for sOrthoSpecies in arg.tOrthoSpecies:
				if sOrthoSpecies not in tListOfRepresentedSpecies:
					bValidClique=False
					try:
						dLost[sOrthoSpecies]+=1
					except KeyError:
						dLost[sOrthoSpecies]=1
			for sDisplaySpecies in arg.tDisplaySpecies:
				if not sDisplaySpecies in dDisplay2Ortholog:
					dDisplay2Ortholog[sDisplaySpecies]={}
				if sDisplaySpecies in tListOfRepresentedSpecies:
					dDisplay2Ortholog[sDisplaySpecies][sCliqueId]="one2one"
				else:
					dDisplay2Ortholog[sDisplaySpecies][sCliqueId]="no ortholog"
		
			if bValidClique:
				tNewListOfGenes=[sRefGeneId]+tListOfGenes
				tNewListOfGenes=sorted(tNewListOfGenes)
				dbListOfGenes=tuple(tNewListOfGenes)
				dClique2Gene[sCliqueId]=dbListOfGenes

	'''
	#Step 6 : OO version
	
	#Pour chaque clique de dCheckInclude : les especes match et ortho sont elles representees ?
	#Si oui, on la met dans un dico dClique2Gene, on peut le laisser comme ca puisque seulement dico simple

	iStepCount=StepIncrease(iStepCount,"Take in count Match/Ortho/Display specificity",0)
	dLost={}
	dClique2Gene={}
	dDisplay2Ortholog={}
	for sRefGeneId in dCheckInclude:
		for sCliqueId in dCheckInclude[sRefGeneId]:
			#tListOfGenes=dCheckInclude[sRefGeneId][sCliqueId]
			tListOfGenes=sCliqueId.otherGenes
			tListOfRepresentedSpecies=[]
			for sSpeciesGeneId in tListOfGenes:
				sSpeciesId=sSpeciesGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
				tListOfRepresentedSpecies.append(sSpeciesId)
		
			bValidClique=True
			for sMatchSpecies in arg.tMatchSpecies:
				if sMatchSpecies not in tListOfRepresentedSpecies:
					bValidClique=False
					try:
						dLost[sMatchSpecies]+=1
					except KeyError:
						dLost[sMatchSpecies]=1
			for sOrthoSpecies in arg.tOrthoSpecies:
				if sOrthoSpecies not in tListOfRepresentedSpecies:
					bValidClique=False
					try:
						dLost[sOrthoSpecies]+=1
					except KeyError:
						dLost[sOrthoSpecies]=1
			for sDisplaySpecies in arg.tDisplaySpecies:
				if not sDisplaySpecies in dDisplay2Ortholog:
					dDisplay2Ortholog[sDisplaySpecies]={}
				if sDisplaySpecies in tListOfRepresentedSpecies:
					dDisplay2Ortholog[sDisplaySpecies][sCliqueId]="one2one"
				else:
					dDisplay2Ortholog[sDisplaySpecies][sCliqueId]="no ortholog"
		
			if bValidClique:
				tNewListOfGenes=[sRefGeneId]+tListOfGenes
				tNewListOfGenes=sorted(tNewListOfGenes)
				dbListOfGenes=tuple(tNewListOfGenes)
				dClique2Gene[sCliqueId]=dbListOfGenes


	'''
	#Step7:retrieve all hits on gene for the jobId
	#NB: if bSingleHitAnalysis is True, retrieve only the best hit
	iStepCount=StepIncrease(iStepCount,"Retrieve all hits for each species",0)
	iCount=1
	sOneTarget = arg.sAbsentValue
	for sSpecies in dGenome2Name.keys():
		if sSpecies== arg.sRefSpecies:
			#Already processed
			continue
		print str(iCount)+"/"+str(len(dGenome2Name)-1)
		iCount+=1
		sSQLcommand="select gene_id,start_match,end_match,sequence_match,strand_match,"+sToolColTarget+",scan_id from "+arg.sVersionDb+"_scan_"+sSpecies+" s where job_id="+arg.sJobId+" "+arg.sregionClause+" and s.match=1;"
		query = SQLquery(sSQLcommand, arg)
		for tRows in query.result:
			sGeneId=tRows[0]
			sSpeciesId=sGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
			dData={"Start":tRows[1],"End":tRows[2],"Sequence":tRows[3],"Strand":tRows[4],sToolColTarget:tRows[5]}
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

			if bAddHit:
				if arg.bSingleHitAnalysis:
					try:
						fStoredScore=float(dGene2Data[sGeneId][dGene2Data[sGeneId].keys()[0]][sToolColTarget])
						if fCurrentScore < fStoredScore:
							dGene2Data[sGeneId]={sScanId:dData}
					except KeyError:
						dGene2Data[sGeneId]={sScanId:dData}
				else:
					try:
						dGene2Data[sGeneId][sScanId]=dData
					except KeyError:
						dGene2Data[sGeneId]={sScanId:dData}
		
			sSpeciesId=sGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId		
			dGenome2GeneSample[sSpeciesId]=sGeneId	
			sOneTarget=tRows[3]

	'''
	#Step 7 : OO version

	#Meme chose qu etape 4 : ajout de chaque hit dans le dico HitDict

	iStepCount=StepIncrease(iStepCount,"Retrieve all hits for each species",0)
	iCount=1
	sOneTarget=arg.sAbsentValue
	for sSpecies in dGenome2Name.keys():
		if sSpecies == arg.sRefSpecies:
			#Already processed
			continue
		print str(iCount)+"/"+str(len(dGenome2Name)-1)
		iCount+=1
		sSQLcommand="select gene_id,start_match,end_match,sequence_match,strand_match,"+sToolColTarget+",scan_id from "+arg.sVersionDb+"_scan_"+sSpecies+" s where job_id="+arg.sJobId+" "+arg.sregionClause+" and s.match=1;"
		#tSQLtable=SubmitMySQLCommand(sSQLcommand,sHost,sDatabase,sUser,sPw)
		query = SQLquery(sSQLcommand,arg)
		for tRows in query.result :
			sGeneId=tRows[0]
			pvalue=float(tRows[5])

			sSpeciesId=sGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
			#dData={"Start":tRows[1],"End":tRows[2],"Sequence":tRows[3],"Strand":tRows[4],sToolColTarget:tRows[5]}
			fCurrentScore=pvalue

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

			if bAddHit:
				sScan_id=tRows[6]
				start=tRows[1]
				end=tRows[2]
				seq=tRows[3]
				strand=tRows[4]	
				hitName=sScan_id+"_"+sGeneId+"_"+start+"_"+end+"_"+seq+"_"+strand+"_"+str(pvalue)
				hitName=Hit(sScan_id,sGeneId,start,end,seq,strand,pvalue)
				if arg.bSingleHitAnalysis:
					try:
						fStoredScore=HitDict[sGeneId][0].pvalue
						if fCurrentScore < fStoredScore:
							HitDict[sGeneId]=[hitName]
					except KeyError:
						HitDict[sGeneId]=[hitName]
				else:
					try:
						HitDict[sGeneId].append(hitName)
					except KeyError:
						HitDict[sGeneId]=[hitName]
		
			sSpeciesId=sGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId		
			dGenome2GeneSample[sSpeciesId]=sGeneId	
			sOneTarget=tRows[3]
	'''
	
	#Step8: Compute and check if hits are homogeneous for MatchSpecies and RefSpecies
	iStepCount=StepIncrease(iStepCount,"Apply homology filter",0)
	tUnalyzedGenes=[]
	if sOneTarget == arg.sAbsentValue:
		print("Warning : no valide results, step skipped")
	elif not arg.bGetResult:
		print("Warning : Count option activated, no homologoy tested, step skipped")
                dClique2Score={}
	elif bCheckHomogeneous:
		dClique2Score={}
	
		iCount=1
		iFinalCount=len(dClique2Gene)
		for sCliqueId in dClique2Gene.keys():
			if iCount%10==0:
				print str(iCount)+"/"+str(iFinalCount)
			iCount+=1
			print sCliqueId
	
			try:
				tTargetedGenes=[X for X in dClique2Gene[sCliqueId] if X[:-11] in tCheckedSpecies]
				dWorkGene2Data={}
				for sGeneId in tTargetedGenes:
					dWorkGene2Data[sGeneId]=dGene2Data[sGeneId]
		
			except KeyError:
				#one or more species have no valid hit for the clique. Clique is removed
				del dClique2Gene[sCliqueId]
				try:
					dLost["Heterogeneous"]+=1
				except KeyError:
					dLost["Heterogeneous"]=1
				continue
			#print dWorkGene2Data,tExceptedPosition,fIdentityMax
	#
                        #Remove Clique if too much combination available
		        iCombination=1
		        bBreak=False
		        for sGeneId in dWorkGene2Data:
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
           		 
			dbTemp=launchHomology4MultiHits(HOMOLOGY_SCRIPT_MODE,dWorkGene2Data,tExceptedPosition,fIdentityMax)
					
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

	'''

	#Step 8 : OO version
	#Pour chaque clique dans dClique2Gene : on ecrit le hit dans un WorkHitDict, qui est la version filtree de HitDict comme dWorkgene2Data est la version filtree de dGene2Data

	iStepCount=StepIncrease(iStepCount,"Apply homology filter",0)
	tUnalyzedGenes=[]
	if sOneTarget==arg.sAbsentValue:
		print("Warning : no valide results, step skipped")
	elif not arg.bGetResult:
		print("Warning : Count option activated, no homology tested, step skipped")
                dClique2Score={}
	elif bCheckHomogeneous:
		dClique2Score={}
	
		iCount=1
		iFinalCount=len(dClique2Gene)
		for sCliqueId in dClique2Gene.keys():
			if iCount%10==0:
				print str(iCount)+"/"+str(iFinalCount)
			iCount+=1
			print sCliqueId
	
			try:
				tTargetedGenes=[X for X in dClique2Gene[sCliqueId] if X[:-11] in tCheckedSpecies]
				dWorkGene2Data={}
				for sGeneId in tTargetedGenes:
					dWorkGene2Data[sGeneId]=dGene2Data[sGeneId]
		
			except KeyError:
				#one or more species have no valid hit for the clique. Clique is removed
				del dClique2Gene[sCliqueId]
				try:
					dLost["Heterogeneous"]+=1
				except KeyError:
					dLost["Heterogeneous"]=1
				continue
			#print dWorkGene2Data,tExceptedPosition,fIdentityMax
	#
                        #Remove Clique if too much combination available
		        iCombination=1
		        bBreak=False
		        for sGeneId in dWorkGene2Data:
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
           		 
			dbTemp=launchHomology4MultiHits(HOMOLOGY_SCRIPT_MODE,dWorkGene2Data,tExceptedPosition,fIdentityMax)
					
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



	#Step9:retriever the ENS id for each species
	iStepCount=StepIncrease(iStepCount,"Retrieve EnsemblId",0)
	dGenome2ENStag={}
	argid=dGenome2GeneSample.values()
	if argid is not None and len(argid)>0:
	
	  sSQLcommand="select stable_id,genome_db_id from "+arg.sVersionDb+"_gene where gene_id="+" or gene_id=".join(dGenome2GeneSample.values())+";"
	  #tSQLtable=SubmitMySQLCommand(sSQLcommand,sHost,sDatabase,sUser,sPw)
	  query = SQLquery(sSQLcommand,arg)
	  for tRows in query.result :
		sENSid=tRows[0]
		sENStag=sENSid[:-11]
		sGenomeId=tRows[1]

		dGenome2ENStag[sGenomeId]=sENStag

	'''
	#Step10:set the global pvalue for each clique

        dClique2UniqPvalue={}
        iStepCount=StepIncrease(iStepCount,"Set the global pvalue for each clique",0)
        for sCliqueId in dClique2Gene:
	  fWorstGlobalPvalue=-1
	  for sGeneId in dClique2Gene[sCliqueId]:
		if sGeneId[:-11] in arg.tCheckedSpecies:
			#Warning: a gene can have no hit, so no occurences in dGene2Data
			try:
				fLocalBestPvalue=1.0
				for sScanId in dGene2Data[sGeneId]:
					fGenePvalue=float(dGene2Data[sGeneId][sScanId][sToolColTarget])
					#print fGenePvalue, fLocalBestPvalue, min(fGenePvalue,fLocalBestPvalue)
					fLocalBestPvalue=min(fGenePvalue,fLocalBestPvalue)
				fWorstGlobalPvalue=max(fLocalBestPvalue,fWorstGlobalPvalue)
			except KeyError:
				continue
	  dClique2UniqPvalue[sCliqueId]=fWorstGlobalPvalue
	'''
	
	#Step10 : OO version

	#Pas bcp de diff si ce n est qu on itere sur HitDict et non sur dGene2Data : on peut acceder direct a la pvalue de la cle en tant qu attribut

	dClique2UniqPvalue={}
        iStepCount=StepIncrease(iStepCount,"Set the global pvalue for each clique",0)
        for sCliqueId in dClique2Gene:
	  fWorstGlobalPvalue=-1
	  for sGeneId in dClique2Gene[sCliqueId]:
		if sGeneId[:-11] in arg.tCheckedSpecies:
			#Warning: a gene can have no hit, so no occurences in dGene2Data
			try:
				fLocalBestPvalue=1.0
				#print HitDict[sGeneId]
				for hit in HitDict[sGeneId] :
					fGenePvalue=hit.pvalue
					#print fGenePvalue, fLocalBestPvalue, min(fGenePvalue,fLocalBestPvalue)
					fLocalBestPvalue=min(fGenePvalue,fLocalBestPvalue)
				fWorstGlobalPvalue=max(fLocalBestPvalue,fWorstGlobalPvalue)
			except KeyError:
				continue
	  dClique2UniqPvalue[sCliqueId]=fWorstGlobalPvalue

	'''
	#Step11:create the final dictionnary
	iStepCount=StepIncrease(iStepCount,"Regroup all results",0)
	dFinal={}
	iCount=1
	for sCliqueId in dClique2Gene:
		if iCount%100==0:
			print str(iCount)+"/"+str(len(dClique2Gene))
		iCount+=1
		#Warning: a clique can have no HGNC annotation
		try:
			dFinal[sCliqueId]={"HGNC":{"full":dClique2HGNC[sCliqueId]["full"],"short":dClique2HGNC[sCliqueId]["short"]},"Content":{}}
		except KeyError:
			dFinal[sCliqueId]={"HGNC":{"full":"no HGNC","short":"no HGNC"},"Content":{}}
		if arg.bCheckHomogeneous:
			dFinal[sCliqueId]["Score"]=dClique2Score[sCliqueId]
		dFinal[sCliqueId]["Pvalue"]=dClique2UniqPvalue[sCliqueId]
		bValidClique=True
		bRefFailure=False
		tMatchFailure=[]
		for sGeneId in dClique2Gene[sCliqueId]:
			#Warning: a gene can have no hit, so no occurences in dGene2Data
			dFinal[sCliqueId]["Content"][sGeneId]={}
			try:
				for sScanId in dGene2Data[sGeneId]:
					dFinal[sCliqueId]["Content"][sGeneId][sScanId]=dGene2Data[sGeneId][sScanId]
			except KeyError:
					sSpeciesId=sGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
					if sSpeciesId == arg.sRefSpecies:
						bValidClique=False
						bRefFailure=True
					elif sSpeciesId in arg.tMatchSpecies:
						bValidClique=False
						tMatchFailure.append(sSpeciesId)
					else:
						dFinal[sCliqueId]["Content"][sGeneId]["noHit"]={"Start":"","End":"","Sequence":"","Strand":"",sToolColTarget:"NA"}					
		if not bValidClique:
			del dFinal[sCliqueId]
			if not bRefFailure:
				for sSpeciesId in tMatchFailure:
					try:
						dLost[sSpeciesId]+=1
					except KeyError:
						dLost[sSpeciesId]=1
	'''

	#Step11 : OO version

	#Meme chose sauf qu on itere sur HitDict et non DGene2Data avec content = attributes ? Possible creation d une classe final_result avec en attribut score, content = Hit

	iStepCount=StepIncrease(iStepCount,"Regroup all results",0)
	dFinal={}
	iCount=1
	for sCliqueId in dClique2Gene:
		if iCount%100==0:
			print str(iCount)+"/"+str(len(dClique2Gene))
		iCount+=1
		#Warning: a clique can have no HGNC annotation
		try:
			dFinal[sCliqueId]={"HGNC":{"full":dClique2HGNC[sCliqueId]["full"],"short":dClique2HGNC[sCliqueId]["short"]},"Content":{}}
		except KeyError:
			dFinal[sCliqueId]={"HGNC":{"full":"no HGNC","short":"no HGNC"},"Content":{}}
		if arg.bCheckHomogeneous:
			dFinal[sCliqueId]["Score"]=dClique2Score[sCliqueId]
		dFinal[sCliqueId]["Pvalue"]=dClique2UniqPvalue[sCliqueId]
		bValidClique=True
		bRefFailure=False
		tMatchFailure=[]
		for sGeneId in dClique2Gene[sCliqueId]:
			#Warning: a gene can have no hit, so no occurences in dGene2Data
			dFinal[sCliqueId]["Content"][sGeneId]={}
			try:
				for hit in HitDict[sGeneId]:
					#dFinal[sCliqueId]["Content"][sGeneId][sScanId]=dGene2Data[sGeneId][sScanId]
					dFinal[sCliqueId]["Content"][sGeneId]=HitDict[sGeneId]
			except KeyError:
					sSpeciesId=sGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
					if sSpeciesId==arg.sRefSpecies:
						bValidClique=False
						bRefFailure=True
					elif sSpeciesId in arg.tMatchSpecies:
						bValidClique=False
						tMatchFailure.append(sSpeciesId)
					else:
						dFinal[sCliqueId]["Content"][sGeneId]["noHit"]={"Start":"","End":"","Sequence":"","Strand":"",sToolColTarget:"NA"}					
		if not bValidClique:
			del dFinal[sCliqueId]
			if not bRefFailure:
				for sSpeciesId in tMatchFailure:
					try:
						dLost[sSpeciesId]+=1
					except KeyError:
						dLost[sSpeciesId]=1


	#Step12:
	iStepCount=StepIncrease(iStepCount,"Compute the ConservationValue",0)
	if arg.bSingleHitAnalysis or arg.bCheckHomogeneous:
		for sCliqueId in dFinal:
			dFinal[sCliqueId]["ConservationValue"]=GetCommonValue(dFinal[sCliqueId]["Content"],tExceptedPosition,tCheckedSpecies)
	else:
		print("Warning : many hits conserved by gene, no ConservationValue computed")
	'''
	#Step13:write the result
	iStepCount=StepIncrease(iStepCount,"Write the results",0)
	if not arg.bGetResult:
		#result = count.json
                if dLost is None:
                    dLost=dict()
                if "Heterogeneous" not in dLost.keys():
                    dLost["Heterogeneous"]=None
		finalRes.dSpeciesCount={"Global":{"cliqueNumber":0,"geneNumber":0,"hitNumber":0,"heterogeneousRejected":dLost["Heterogeneous"]},"Reference":{"geneNumber":0,"hitNumber":0},"MatchSpecies":{},"OrthoSpecies":{},"DisplaySpecies":{}}
		dSpeciesClassification={}
		dSpeciesClassification[arg.sRefSpecies]={"type":"Reference","name":dGenome2Name[arg.sRefSpecies]}
		for sSpeciesId in arg.tMatchSpecies:
			dSpeciesClassification[sSpeciesId]={"type":"MatchSpecies","name":dGenome2Name[sSpeciesId]}
			try:
				iLostNumber=dLost[sSpeciesId]
			except KeyError:
				iLostNumber=0
			finalRes.dSpeciesCount["MatchSpecies"][dGenome2Name[sSpeciesId]]={"geneNumber":0,"hitNumber":0,"lostClique":iLostNumber}
		for sSpeciesId in arg.tOrthoSpecies:
			dSpeciesClassification[sSpeciesId]={"type":"OrthoSpecies","name":dGenome2Name[arg.sRefSpecies]}
			try:
				iLostNumber=dLost[sSpeciesId]
			except KeyError:
				iLostNumber=0
			finalRes.dSpeciesCount["OrthoSpecies"][dGenome2Name[sSpeciesId]]={"geneNumber":0,"hitNumber":0,"lostClique":iLostNumber}
		for sSpeciesId in arg.tDisplaySpecies:
			dSpeciesClassification[sSpeciesId]={"type":"DisplaySpecies","name":dGenome2Name[arg.sRefSpecies]}
			finalRes.dSpeciesCount["DisplaySpecies"][dGenome2Name[sSpeciesId]]={"geneNumber":0,"hitNumber":0}

		for sCliqueId in dFinal:
			finalRes.dSpeciesCount["Global"]["cliqueNumber"]+=1
			for sGeneId in dFinal[sCliqueId]["Content"]:

				finalRes.dSpeciesCount["Global"]["geneNumber"]+=1
				sSpeciesId=sGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
				sSpeciesName=dGenome2Name[sSpeciesId]
				iHitNumber=len(dFinal[sCliqueId]["Content"][sGeneId])
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
	else:
		#Add BestPvalue column
		bBestPvalue=False
		if not bSingleHitAnalysis and not bCheckHomogeneous:
			bBestPvalue=True
	
		#Create headline of the file
		sHeaderLine="HGNC long\tHGNC short\tBiotype\t"
		if len(tExceptedPosition)==0:
			sExceptedTag=""
		else:
			sExceptedTag="[excepted-"+"-".join(tExceptedPosition)+"]"
		if bCheckHomogeneous:
			sHeaderLine+="HomologyScore"+sExceptedTag+"\t"
		if bCheckHomogeneous or bSingleHitAnalysis:
			sHeaderLine+="ConservationValue"+sExceptedTag+"\t"
		sHeaderLine+="WorstPvalue\t"
		for sSpeciesId in tAllSpecies:
			sSpeciesShortId=""
			if sSpeciesId in dGenome2ENStag:
				sSpeciesShortId=dGenome2ENStag[sSpeciesId][:-1]
				sSpeciesShortId=sSpeciesShortId.replace("ENS","")
			if sSpeciesShortId=="":
				sSpeciesShortId="HSA"
		
			sSpeciesCategory=sSpeciesShortId+"."
			if sSpeciesId==sRefSpecies:
				sSpeciesCategory+="Ref"
			elif sSpeciesId in tMatchSpecies:
				sSpeciesCategory+="Match"
			elif sSpeciesId in tOrthoSpecies:
				sSpeciesCategory+="Ortho"
			elif sSpeciesId in tDisplaySpecies:
				sSpeciesCategory+="Display"
			else:
				exit("Unknow speciesId "+sSpeciesId+", process broken")
			
			sHeaderLine+=sSpeciesCategory
			if bBestPvalue:
				sHeaderLine+="\t"+sSpeciesShortId+".BestPvalue"
			sHeaderLine+="\t"+sSpeciesShortId+".gene\t"+sSpeciesShortId+"."+sToolColTarget+"\t"+sSpeciesShortId+".coord\t"+sSpeciesShortId+".sequence\t"
			if sSpeciesId in tDisplaySpecies:
				sHeaderLine+=sSpeciesShortId+".ortholog\t"

		FILE=open(sResultFileName+".tsv","w")
		FILE.write(sHeaderLine+"\n")	

		#Write each line, 1 by valide clique
		iCount=1
		iWritedClique=0
		for sCliqueId in dFinal:
			if iCount%100==0:
				print str(iCount)+"/"+str(len(dFinal))
			iCount+=1

			#HGNC annotation
			sCurrentLine=dFinal[sCliqueId]["HGNC"]["full"]+"\t"+dFinal[sCliqueId]["HGNC"]["short"]+"\t"
			#Biotype annotation
			sCurrentLine+=dClique2Biotype[sCliqueId]+"\t"
			#Identity Score
			if bCheckHomogeneous:
				sCurrentLine+=str(round(dFinal[sCliqueId]["Score"],2))+"\t"
			if bCheckHomogeneous or bSingleHitAnalysis:
				sCurrentLine+=str(dFinal[sCliqueId]["ConservationValue"])+"\t"
			#Pvalue Score
			sCurrentLine+=str(dFinal[sCliqueId]["Pvalue"])+"\t"

			#For each species
			for sSpeciesId in tAllSpecies:
				asens=0
				if sSpeciesId in dGenome2ENStag:
					asens=1	
				else:
					asens=0
					dGenome2ENStag[sSpeciesId]=""
				#CategorySpacer
				sCurrentLine+="\t"
				bNoGeneValide=True
				#Search the affiliate gene
				for sGeneId in dFinal[sCliqueId]["Content"]:
					if sGeneId[:-11]!=sSpeciesId: #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
						continue

					bNoGeneValide=False
					#two case : have at least one valide hit or have no hit
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

	'''
	#step13 : OO version
	iStepCount=StepIncrease(iStepCount,"Write the results",0)
	if not arg.bGetResult:
		#result = count.json
                if dLost is None:
                    dLost=dict()
                if "Heterogeneous" not in dLost.keys():
                    dLost["Heterogeneous"]=None
		finalRes.dSpeciesCount={"Global":{"cliqueNumber":0,"geneNumber":0,"hitNumber":0,"heterogeneousRejected":dLost["Heterogeneous"]},"Reference":{"geneNumber":0,"hitNumber":0},"MatchSpecies":{},"OrthoSpecies":{},"DisplaySpecies":{}}
		dSpeciesClassification={}
		dSpeciesClassification[arg.sRefSpecies]={"type":"Reference","name":dGenome2Name[arg.sRefSpecies]}
		for sSpeciesId in arg.tMatchSpecies:
			dSpeciesClassification[sSpeciesId]={"type":"MatchSpecies","name":dGenome2Name[sSpeciesId]}
			try:
				iLostNumber=dLost[sSpeciesId]
			except KeyError:
				iLostNumber=0
			finalRes.dSpeciesCount["MatchSpecies"][dGenome2Name[sSpeciesId]]={"geneNumber":0,"hitNumber":0,"lostClique":iLostNumber}
		for sSpeciesId in arg.tOrthoSpecies:
			dSpeciesClassification[sSpeciesId]={"type":"OrthoSpecies","name":dGenome2Name[arg.sRefSpecies]}
			try:
				iLostNumber=dLost[sSpeciesId]
			except KeyError:
				iLostNumber=0
			finalRes.dSpeciesCount["OrthoSpecies"][dGenome2Name[sSpeciesId]]={"geneNumber":0,"hitNumber":0,"lostClique":iLostNumber}
		for sSpeciesId in arg.tDisplaySpecies:
			dSpeciesClassification[sSpeciesId]={"type":"DisplaySpecies","name":dGenome2Name[arg.sRefSpecies]}
			finalRes.dSpeciesCount["DisplaySpecies"][dGenome2Name[sSpeciesId]]={"geneNumber":0,"hitNumber":0}

		for sCliqueId in dFinal:
			finalRes.dSpeciesCount["Global"]["cliqueNumber"]+=1
			for sGeneId in dFinal[sCliqueId]["Content"]:

				finalRes.dSpeciesCount["Global"]["geneNumber"]+=1
				sSpeciesId=sGeneId[:-11] #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
				sSpeciesName=dGenome2Name[sSpeciesId]
				iHitNumber=len(dFinal[sCliqueId]["Content"][sGeneId])
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
	else:
		#Add BestPvalue column
		bBestPvalue=False
		if not bSingleHitAnalysis and not bCheckHomogeneous:
			bBestPvalue=True
	
		#Create headline of the file
		sHeaderLine="HGNC long\tHGNC short\tBiotype\t"
		if len(tExceptedPosition)==0:
			sExceptedTag=""
		else:
			sExceptedTag="[excepted-"+"-".join(tExceptedPosition)+"]"
		if bCheckHomogeneous:
			sHeaderLine+="HomologyScore"+sExceptedTag+"\t"
		if bCheckHomogeneous or bSingleHitAnalysis:
			sHeaderLine+="ConservationValue"+sExceptedTag+"\t"
		sHeaderLine+="WorstPvalue\t"
		for sSpeciesId in tAllSpecies:
			sSpeciesShortId=""
			if sSpeciesId in dGenome2ENStag:
				sSpeciesShortId=dGenome2ENStag[sSpeciesId][:-1]
				sSpeciesShortId=sSpeciesShortId.replace("ENS","")
			if sSpeciesShortId=="":
				sSpeciesShortId="HSA"
		
			sSpeciesCategory=sSpeciesShortId+"."
			if sSpeciesId==sRefSpecies:
				sSpeciesCategory+="Ref"
			elif sSpeciesId in tMatchSpecies:
				sSpeciesCategory+="Match"
			elif sSpeciesId in tOrthoSpecies:
				sSpeciesCategory+="Ortho"
			elif sSpeciesId in tDisplaySpecies:
				sSpeciesCategory+="Display"
			else:
				exit("Unknow speciesId "+sSpeciesId+", process broken")
			
			sHeaderLine+=sSpeciesCategory
			if bBestPvalue:
				sHeaderLine+="\t"+sSpeciesShortId+".BestPvalue"
			sHeaderLine+="\t"+sSpeciesShortId+".gene\t"+sSpeciesShortId+"."+sToolColTarget+"\t"+sSpeciesShortId+".coord\t"+sSpeciesShortId+".sequence\t"
			if sSpeciesId in tDisplaySpecies:
				sHeaderLine+=sSpeciesShortId+".ortholog\t"

		FILE=open(sResultFileName+".tsv","w")
		FILE.write(sHeaderLine+"\n")	

		#Write each line, 1 by valide clique
		iCount=1
		iWritedClique=0
		for sCliqueId in dFinal:
			if iCount%100==0:
				print str(iCount)+"/"+str(len(dFinal))
			iCount+=1

			#HGNC annotation
			sCurrentLine=dFinal[sCliqueId]["HGNC"]["full"]+"\t"+dFinal[sCliqueId]["HGNC"]["short"]+"\t"
			#Biotype annotation
			sCurrentLine+=dClique2Biotype[sCliqueId]+"\t"
			#Identity Score
			if bCheckHomogeneous:
				sCurrentLine+=str(round(dFinal[sCliqueId]["Score"],2))+"\t"
			if bCheckHomogeneous or bSingleHitAnalysis:
				sCurrentLine+=str(dFinal[sCliqueId]["ConservationValue"])+"\t"
			#Pvalue Score
			sCurrentLine+=str(dFinal[sCliqueId]["Pvalue"])+"\t"

			#For each species
			for sSpeciesId in tAllSpecies:
				asens=0
				if sSpeciesId in dGenome2ENStag:
					asens=1	
				else:
					asens=0
					dGenome2ENStag[sSpeciesId]=""
				#CategorySpacer
				sCurrentLine+="\t"
				bNoGeneValide=True
				#Search the affiliate gene
				for sGeneId in dFinal[sCliqueId]["Content"]:
					if sGeneId[:-11]!=sSpeciesId: #In geneId, the last 11 numbers are the numerical part of ENS id. The part before is the speciesId
						continue

					bNoGeneValide=False
					#two case : have at least one valide hit or have no hit
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

		#Iterer sur le dico result version classe
	'''
mainprocess()

########fix :cleaning temp data
cmd=["bash","-c", "rm -f ./TempCompa*" ]
ExecuteBashCommand(cmd)
cmd=["bash","-c", "rm -f ./datatmp*" ]
ExecuteBashCommand(cmd)
