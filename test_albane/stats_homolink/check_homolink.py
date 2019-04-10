with open("global_homolink.txt","r") as homolink :
	homolink.readline()
	ct = 0
	for line in homolink :
		ct += 1
		species_id=line.split("\t")
		id1=species_id[0][0:3]
		id2=species_id[1][0:3]
		if id1 == id2 :
			print(line)
		if ct%1000000 == 0 :
			print(ct)

