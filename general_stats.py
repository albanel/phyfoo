import matplotlib.pyplot as plt
import time
import pymysql
import sys

class study_file() :
	def __init__(self,file) :
		self.start_time = time.time()
		self.file = open(file,"r")
		self.filename = file
		self.outpout_name = file.replace(".txt","_statistics.txt")
		self.output = open(self.outpout_name,"w")
	
	def get_dico_specie_to_id(self) :
		conn = pymysql.connect(host='192.168.1.198',  user='orthocis3', passwd='r77ujtt9', db='orthocis3')
		cur = conn.cursor()
		self.dico_id_to_speciename = {}
		for specie_id in self.specie_id_to_link_number.keys() :
			req = "SELECT name FROM current_genome_db WHERE genome_db_id = "+specie_id+" ;"
			cur.execute(req)
			for elt in cur :
				self.dico_id_to_speciename[specie_id] = elt[0]
				if elt[0] == "" :
					print("No name for specie "+specie_id)
		cur.close()
		conn.close()
		#print(self.dico_id_to_speciename)

	def show_attributes(self) :
		# General info : number of genes and links in the network
		self.output.write("The network has "+str(self.number_of_genes)+" genes and "+str(self.link_number)+" links between them."+"\n")
		self.output.write("The network involves "+str(self.species_number)+ " species."+"\n")
		cumuled_gene_number = 0
		cumuled_link_number = 0
		# Species info : number of gene and links for each specie
		for specie in self.specie_id_to_gene_number.keys() :
			self.output.write("The specie "+str(self.dico_id_to_speciename[specie])+" has "+str(self.specie_id_to_gene_number[specie])+" total genes and "+str(self.specie_id_to_link_number[specie])+" total links in the network for an average connectivity of "+str(self.specie_connexion[specie])+"\n")
			cumuled_gene_number += self.specie_id_to_gene_number[specie]
			cumuled_link_number += self.specie_id_to_link_number[specie]
		# Species with the three most import gene and links number
		ordered_gene_numbers = sorted(self.specie_id_to_gene_number.values(), reverse = True) 
		ordered_link_numbers = sorted(self.specie_id_to_link_number.values(), reverse = True)
		ordered_gene_numbers2 = sorted(self.specie_id_to_gene_number.values(), reverse = False)
		ordered_link_numbers2 = sorted(self.specie_id_to_link_number.values(), reverse = False)
		self.output.write("Species with the most important gene number : "+"\n")
		for specie in self.specie_id_to_gene_number.keys() :
			if self.specie_id_to_gene_number[specie] in ordered_gene_numbers[:3] :
				self.output.write("\t"+str(self.dico_id_to_speciename[specie])+", gene number : "+str(self.specie_id_to_gene_number[specie])+"\n")
		self.output.write("Species with the most important link number : "+"\n")
		for specie in self.specie_id_to_link_number.keys() :
			if self.specie_id_to_link_number[specie] in ordered_link_numbers[:3] :
				self.output.write("\t"+str(self.dico_id_to_speciename[specie])+", link number : "+str(self.specie_id_to_link_number[specie])+"\n")

		self.output.write("Species with the least important gene number : "+"\n")
		for specie in self.specie_id_to_gene_number.keys() :
			if self.specie_id_to_gene_number[specie] in ordered_gene_numbers2[:3] :
				self.output.write("\t"+str(self.dico_id_to_speciename[specie])+", gene number : "+str(self.specie_id_to_gene_number[specie])+"\n")
		self.output.write("Species with the least important link number : "+"\n")
		for specie in self.specie_id_to_link_number.keys() :
			if self.specie_id_to_link_number[specie] in ordered_link_numbers2[:3] :
				self.output.write("\t"+str(self.dico_id_to_speciename[specie])+", link number : "+str(self.specie_id_to_link_number[specie])+"\n")
		# Average gene and link number
		avg_gene_number = cumuled_gene_number/self.species_number 
		avg_link_number = cumuled_link_number/self.species_number
		self.output.write("Average gene number per specie : "+str(avg_gene_number)+"\n")
		self.output.write("Average link number per specie : "+str(avg_link_number)+"\n")
		# The most connected genes are...
		self.output.write("Statistics for all genes in the network :"+"\n")
		cumuled_connectivity = sum(self.gene_connexion_distribution)
		avg_connectivity = cumuled_connectivity/self.number_of_genes
		self.output.write("\t"+"The average connectivity for each gene in the network is "+str(avg_connectivity)+"\n")
		self.output.write("\t"+"The max connection in the gene network is "+str(max(self.gene_connexion_distribution))+". This value is found for "+str(self.gene_connexion_distribution.count(max(self.gene_connexion_distribution)))+" genes ("+str(self.gene_connexion_distribution.count(max(self.gene_connexion_distribution))/self.number_of_genes)+"% of genes)."+"\n") 		   
		self.output.close()
		print("The statistics where successfully recorded in the file "+self.outpout_name+".")
		self.end_time = time.time()
		print("Time ellapsed : "+str(self.end_time-self.start_time)+" seconds.")

	def get_histogram(self) :
		gene_connexion_distribution = [] 
		for value in self.gene_links.values() :
			gene_connexion_distribution.append(len(value))
		self.gene_connexion_distribution = gene_connexion_distribution
		#print(self.gene_connexion_distribution)
		plt.hist(self.gene_connexion_distribution, bins=30)
		# naming the x axis 
		plt.ylabel("Nombre de gènes")
		# naming the y axis 
		plt.xlabel("Connectivité")
		# giving a title to my graph 
		plt.title("Connectivité des gènes dans la table global_homolink")
                # save the plot
		title_save = self.filename.replace(".txt","")+"_histogram.png"
		plt.savefig(title_save)
		print("The histogram was successfully saved under the name "+title_save+".")		

	def all_analysis(self) :
		print("Statistical analysis : start")
		self.file.readline()
		self.gene_links = {}
		self.specie_id_to_gene_number = {}
		self.specie_id_to_link_number = {}
		self.gene_connexion = {}
		self.specie_connexion = {}
		self.link_number = 0
		for line in self.file :
			self.link_number += 1
			line = line.replace("\n","")
			line = line.split("\t") ## Linked genes are opposite in each line of the file, separated by a tab
			gene1 = line[0]
			gene2 = line[1]
			if gene1 in self.gene_links.keys() : ## If the gene is already in the dictionnary, the opposite is added
				self.gene_links[gene1].append(gene2)
			else :
				self.gene_links[gene1]=[gene2] ## Or its initialized with a list with the linked gene inside
			if gene2 in self.gene_links.keys() : ## Same for the opposite gene
				self.gene_links[gene2].append(gene1)
			else :
				self.gene_links[gene2]=[gene1]
		self.number_of_genes = len(self.gene_links.keys())
		print("Calculating the number of links per gene...")
		for gene in self.gene_links.keys() :
			if gene in self.gene_connexion.keys() :
				self.gene_connexion[gene] += len(self.gene_links[gene])
			else :
				self.gene_connexion[gene] = len(self.gene_links[gene])
		print("Calculating the number of links per specie...")
		for gene in self.gene_links.keys() :
			specie_id = gene[0:3]
			if specie_id in self.specie_id_to_link_number.keys() :
				self.specie_id_to_link_number[specie_id] += len(self.gene_links[gene])
			else :
				self.specie_id_to_link_number[specie_id] = len(self.gene_links[gene])
		print("Calculating the number of genes per specie...")
		for gene in self.gene_links.keys() :
			specie_id = gene[0:3]
			if specie_id in self.specie_id_to_gene_number.keys() :
				self.specie_id_to_gene_number[specie_id] += 1 
			else :
				self.specie_id_to_gene_number[specie_id] = 1
		print("Calculating the connectivity for each specie...")
		for specie in self.specie_id_to_gene_number.keys() :
			self.specie_connexion[specie] = self.specie_id_to_link_number[specie]/self.specie_id_to_gene_number[specie]
		self.species_number = len(self.specie_id_to_link_number.keys())
		print("Plotting the histogram...")
		self.get_histogram()
		self.file.close()
		self.get_dico_specie_to_id()
		self.show_attributes()
		print("Statistical analysis : end")		
		
homolink = study_file(sys.argv[1])

homolink.all_analysis()

#global_homolink = study_file("global_homolink.txt")

#global_homolink.all_analysis()

#To do

#Number of links, stats max and min number of links, which gene, average

#Stats by specie : max and min number of links, which specie, average

#Stats about link type
