import pymysql

conn = pymysql.connect(host='192.168.1.198',  user='orthocis3', passwd='r77ujtt9', db='orthocis3')
cur = conn.cursor()

dict = eval(open("dSpeciesToLink_global.txt").read())

dict2={}

for genome_id in dict.keys() :
	req = "SELECT name FROM current_genome_db WHERE genome_db_id = "+genome_id+" ;"
	cur.execute(req)
	for specie in cur :
		if specie[0] == "" :
			print("Pas de non pour l'espèce "+genome_id)
		specie=specie[0]
		value=dict[genome_id]
		dict2[specie]=value

#print(dict2)

fichier_sortie=open("species_name_to_link_global.txt","w")

for valeur in sorted(dict2.values()) :
	for genome_id in dict2.keys() :
		if dict2[genome_id] == valeur :
			fichier_sortie.write(str(genome_id)+" "+str(valeur)+"\n")
fichier_sortie.close()

cur.close()
conn.close()

