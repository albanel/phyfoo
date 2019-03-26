with open("global_homolink.txt","r") as homolink :
	homolink.readline()
	for line in homolink :
		species_id=line.split("\t")
		id1=species_id[0][0:3]
		id2=species_id[1][0:3]
		if id1 == id2 :
			print(line)

