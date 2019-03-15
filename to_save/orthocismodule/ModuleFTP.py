from ftplib import FTP
import re
import sys
import os

def download(filename,oudir,ftp,newName="ThisIsADefaultValueIfYouWillNotRenameTheDowloadedFile"):
	if newName=="ThisIsADefaultValueIfYouWillNotRenameTheDowloadedFile":
		newName=filename
	# Open the file for writing in binary mode
	print 'Opening local file ' + filename
	file = open(oudir+"/"+newName, 'wb')

	print 'Getting ' + filename
	ftp.retrbinary('RETR %s' % filename, file.write)

	# Clean up time
	print 'Closing file ' + filename
	file.close()

def listContent(cdir,ftp):
	ftp.cwd(cdir)
	lst=[];
	for name in ftp.nlst():
		lst.append(name)
	return lst


'''
def filterContent(lst,exp):
	filt=[]
	p = re.compile(exp)
	print p
	for name in lst:
		#print "listing: " + name
		if p.match(name):
			filt.append(name)
	return filt
'''
##albane on peut telecharger les especes voulues, et pas seulement une ou toutes
def filterContent(lst,exp):
	#print exp
        filt=[]
	liste_espece=exp.split(",")
	for espece in liste_espece :
		#print espece
		p = re.compile(espece)
        	for name in lst:
                	#print "listing: " + name
                	if p.match(name):
                        	filt.append(name)
        return filt


def updir(ftp):
	ftp.cwd("../")
 


def ftpSpeciesDataDownload(doDownload,ftpurl,oudir,localfastadir,speciesfilter,version):
	log=[]
	ftp = FTP(ftpurl) 
	ftp.login()  
	 
	print "---------------------"
	print "-Download Core"
	rootd="/pub"
	rcont0=listContent(rootd,ftp)
	dircontent0=filterContent(rcont0,'release-'+version)
	for f0 in dircontent0:
		print "release-X"
		rcont=listContent(f0,ftp)
		dircontent=filterContent(rcont,'mysql')
		for f1 in dircontent:
			print "mysql"
			dircontent2=filterContent(listContent(f1,ftp),'.*core.*')
			dircontent2=filterContent(dircontent2,speciesfilter)
			for f2 in dircontent2:
				print "core"
			 	#print "\tFOUND2: " + f2
				dircontent3=filterContent(listContent(f2,ftp),'transcript.txt.gz|assembly.txt.gz|coord_system.txt.gz|seq_region.txt.gz|gene.txt.gz')
				for f3 in dircontent3:
					fil=rootd+"/"+f1+"/"+f2+"/"+f3
			 		print "\t\tFOUND: " +fil

	#				print(oudir+"/"+f2+"/"+f3)
					currentOutdir=oudir+"/"+f2
					try:
						os.makedirs(currentOutdir)
					except OSError :
						print(currentOutdir+" already exist")
				
					if(doDownload):
						download(f3,currentOutdir,ftp)
						log.append(fil)
				updir(ftp)
			updir(ftp)
		updir(ftp)
		ftp.quit() 

	ftp = FTP(ftpurl) 
	ftp.login()  
	print "---------------------"
	print "-Download Fasta"
	rootd="/pub"
	rcont0=listContent(rootd,ftp)
	dircontent0=filterContent(rcont0,'release-'+version)
	for f0 in dircontent0:
		print "release-X"
		rcont=listContent(f0,ftp)
		dircontent=filterContent(rcont,'fasta')
		for f1 in dircontent:
			print "fasta"
			dircontent2=filterContent(listContent(f1,ftp),'.*')
			dircontent2=filterContent(dircontent2,speciesfilter)
			for f2 in dircontent2:
				dircontent3=filterContent(listContent(f2,ftp),'dna')
				for f3 in dircontent3:
			 		#print "\tFOUND3: " + f3
					dircontent4=filterContent(listContent(f3,ftp),'.*dna.toplevel.fa.gz|.*.dna.primary_assembly.fa.gz')
					for f4 in dircontent4:
						fil=rootd+"/"+f1+"/"+f2+"/"+f3+"/"+f4
			 			print "\t\tFOUND: " +fil
						if(doDownload):
							download(f4,localfastadir,ftp)
							log.append(fil)
					updir(ftp)
				updir(ftp)
			updir(ftp)
		updir(ftp)
	 	ftp.quit() 
	 	
	return log


def ftpSpeciesRelationDownload(doDownload,ftpurl,oudir,speciesfilter,version):
	log=[]
	print ftpurl
	ftp = FTP(ftpurl) 
	ftp.login()  
	 
	print "---------------------"
	rootd="/pub"
	rcont0=listContent(rootd,ftp)
	print rootd
	dircontent0=filterContent(rcont0,'release-'+version)
	for f0 in dircontent0:
		print f0
		rcont=listContent(f0,ftp)
		dircontent=filterContent(rcont,'mysql')
		for f1 in dircontent:
			 print "FOUND1: " + f1
			 
			 dircontent2=filterContent(listContent(f1,ftp),'ensembl_compara_.*')
			 for f2 in dircontent2:
				#version=f2.split("_")[-1]
				dircontent3=filterContent(listContent(f2,ftp),'genome_db.txt.gz|homology.txt.gz|homology_member.txt.gz|gene_member.txt.gz|method_link_species_set.txt.gz|species_set.txt.gz|method_link.txt.gz')

				for f3 in dircontent3:
					fil=rootd+"/"+f1+"/"+f2+"/"+f3
			 		print "\t\tFOUND: " +fil
					currentOutdir=oudir+"/"+f2
					try:
						os.makedirs(currentOutdir)
					except OSError :
						print(currentOutdir+" already exist")
				
					if(doDownload):
						download(f3,currentOutdir,ftp,version+"_"+f3)
						log.append(fil)
				updir(ftp)
			 updir(ftp)
		updir(ftp)
	updir(ftp)
	ftp.quit()
	return log

def ftpHGNCfileDownload(doDownload,ftpurl,currentOutdir,pathFile):
	log=[]
	ftp=FTP(ftpurl)
	ftp.login()
	
	
	try:
		os.makedirs(currentOutdir)
	except OSError :
		print(currentOutdir+" already exist")	
	
	if doDownload:
		download(pathFile,currentOutdir,ftp,"HGNC_complete_set.txt.gz")
		log.append(pathFile)	
	ftp.quit()
	return log






