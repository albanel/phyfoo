import pymysql

conn = pymysql.connect(host='192.168.1.198',  user='orthocis3', passwd='r77ujtt9', db='orthocis3')
cur = conn.cursor()

tables_to_drop=["current_assembly","current_coord","current_dna","current_gene_","current_scan","seq_region","current_sequence","current_transcript"]

tables_to_drop2=[]

for table in tables_to_drop :
	req1 = "SHOW TABLES LIKE " + "'%" + table + "%'" + " ;"
	print(req1)
	cur.execute(req1)
	for elt in cur : 
		tables_to_drop2.append(elt)



for elt in tables_to_drop2 :
   elt = str(elt)
   elt = elt.replace("(","")
   elt = elt.replace(")","")
   elt = elt.replace(",","")
   elt = elt.replace("'","")
   req = "DROP TABLE " + elt +" ;"
   print(req)
   cur.execute(req)

cur.close()

conn.close()

