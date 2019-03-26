with open("global_homolink.txt","r") as homolink :
	homolink.readline()
	gene_links={}
	#dSpeciesToLink={}
	for line in homolink :
		line = line.replace("\n","")
		line = line.split("\t")
		gene1 = line[0]
		gene2 = line[1]
		if gene1 in gene_links.keys() :
			gene_links[gene1].append(gene2)
		else :
			gene_links[gene1]=[]
			gene_links[gene1].append(gene2)
		if gene2 in gene_links.keys() :
			gene_links[gene2].append(gene1)
		else :
			gene_links[gene2]=[]
			gene_links[gene2].append(gene1)
	link_length=[]
	for gene in gene_links.keys() : 
		link_length.append(len(gene_links[gene]))
	sorted_list=sorted(link_length,reverse=True)
	print(max(sorted_list))
	print(sorted_list[0:100])
	for gene in gene_links.keys() :
		if len(gene_links[gene]) == max(sorted_list) :
			print("Le gène "+gene+" a le plus de connexion avec ses voisins dans le réseau global.")
	#for gene in gene_links.keys() :
	#	specie=gene[0:3]
	#	if specie in dSpeciesToLink.keys() :
	#		dSpeciesToLink[specie] += len(gene_links[gene])
	#	else :
	#		dSpeciesToLink[specie] = 0
	#		dSpeciesToLink[specie] += len(gene_links[gene])

#with open("dSpeciesToLink.txt","w") as dSpeciesToLinkFile :		
#	dSpeciesToLinkFile.write(str(dSpeciesToLink))

'''
#homolink=open("homolink_distribution.txt","w")
max_homolink=0
nbre_genes=0
nbre_liens_total=0
#dico_nbre_liens={}
for gene in gene_links.keys() :
	nbre_genes+=1
	nbre_liens = len(gene_links[gene])
	nbre_liens_total+=nbre_liens
	if nbre_liens > max_homolink :
		max_homolink = nbre_liens

#	homolink.write(str(nbre_liens)+"\t")
#        if nbre_liens in dico_nbre_liens.keys() :
#                dico_nbre_liens[nbre_liens]+=1
#        else :
#                dico_nbre_liens[nbre_liens]=1

#homolink.close()

print("Le nombre de liens au max dans homolink est de "+str(max_homolink))
print("Le nombre moyen de lien par gène est de "+str(nbre_liens_total/nbre_genes))

with open("global_homolink.txt","r") as homolink :
	homolink.readline()
	gene_links={}
	ct = 0
	for line in homolink :
		line = line.replace("\n","") 
		line = line.split("\t")
		gene1 = line[0]
		gene2 = line[1]
		if gene1 in gene_links.keys() :
			gene_links[gene1].append(gene2)
		else :
			gene_links[gene1]=[]
			gene_links[gene1].append(gene2)
		if gene2 in gene_links.keys() :
			gene_links[gene2].append(gene1)
		else :
			gene_links[gene2]=[]
			gene_links[gene2].append(gene1)

	for gene in gene_links.keys() :
		specie=gene[0:3]
		if specie in dSpeciesToLink.keys() :
			dSpeciesToLink[specie] += len(gene_links[gene])
		else :
			dSpeciesToLink[specie] = 0
			dSpeciesToLink[specie] += len(gene_links[gene])

with open("dSpeciesToLink_global.txt","w") as dSpeciesToLinkFile :
        dSpeciesToLinkFile.write(str(dSpeciesToLink))


#homolink=open("global_homolink_distribution.txt","w")
max_global_homolink=0
nbre_global_genes=0
nbre_liens_global_total=0
#dico_nbre_liens={}
for gene in gene_links.keys() :
	nbre_global_genes+=1
	nbre_liens=len(gene_links[gene])
	nbre_liens_global_total+=nbre_liens
	if nbre_liens > max_global_homolink :
		max_global_homolink = nbre_liens

#        homolink.write(str(nbre_liens)+"\t")
#        if nbre_liens in dico_nbre_liens.keys() :
#                dico_nbre_liens[nbre_liens]+=1
#        else :
#                dico_nbre_liens[nbre_liens]=1

#homolink.close()
print("Le nombre de liens au max dans global_homolink est de "+str(max_global_homolink))
print("Le nombre moyen de lien par gène est de "+str(nbre_liens_global_total/nbre_global_genes))
'''
