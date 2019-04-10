import pymysql
import sys
import os


# This programm enables to know the type of homology between two genes whose id (like in homolink tables) are given in arguments

def get_description(gene1,gene2) :

	conn = pymysql.connect(host='192.168.1.198',  user='orthocis3', passwd='r77ujtt9', db='orthocis3')
	cur = conn.cursor()

	# First of all we get the stable id from the transformed gene_id thanks to the table current_gene

	req1 = "select a.stable_id from current_gene a, current_member b where global_id="+gene1+" or global_id="+gene2+" limit 2 ;"

	#print(req1)
	cur.execute(req1)
	stable_ids = []
	for elt in cur : 
		stable_ids.append(elt[0])

	#print(stable_ids)

	member_ids = []

	# For each stable id, we want to get the member_id in the table member

	for stable_id in stable_ids :
		req2 = "select member_id from current_member where stable_id = '"+stable_id+"' ;"
		#print(req2)
		cur.execute(req2)
		for elt in cur :
			member_ids.append(elt[0])

	#print(member_ids)

	homolist = []

	# For each member_id, we want to het the homology_id in the table homology_member

	for member_id in member_ids :
		req3 = "select homology_id from current_homology_member where member_id = '"+str(member_id)+"' ;"
		cur.execute(req3)	
		for elt in cur :
			homolist.append(elt[0])

	# We have the homolist with all the homology_id for both genes. The value with links both is the only value which appears twice

	for elt in homolist :
		if homolist.count(elt) == 2 :
			#print("2 genes united under the homology id "+str(elt))
			req4 = "select description from current_homology where homology_id = '"+str(elt)+"' ;"
			cur.execute(req4)
			for elt2 in cur :
				#print("Type of this homology : "+str(elt2[0]))
				return(str(elt2[0])+"\n")
			break
	cur.close()
	conn.close()


different_types_of_links = []

filename = sys.argv[1]

with open(filename,"r") as homolink :
	homolink.readline()
	for line in homolink :
		line = line.replace("\n","")
		line = line.split("\t")
		gene1 = line[0]
		gene2 = line[1]
		if get_description(gene1,gene2) not in different_types_of_links :
			different_types_of_links.append(get_description(gene1,gene2))
			os.system("echo 'new type of link in "+filename+" : "+get_description(gene1,gene2)+"' | mail albane.lysiak@etudiant.univ-rennes1.fr")
	
os.system("echo 'Finished for "+filename+ "' | mail albane.lysiak@etudiant.univ-rennes1.fr")







