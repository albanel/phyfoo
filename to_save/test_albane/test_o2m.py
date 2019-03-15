fichier=open("test_o2m.txt","r")

for ligne in fichier :
	print ligne
	sep=ligne.split("\t")
	print sep[0]
	print sep[1]
	print sep[2]
	print sep[3:]
	break

fichier.close()
