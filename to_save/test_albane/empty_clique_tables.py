import pymysql

conn = pymysql.connect(host='192.168.1.198',  user='orthocis3', passwd='r77ujtt9', db='orthocis3')
cur = conn.cursor()

tables_to_empty=["clique","compo"]

tables_to_empty2=[]

for pattern in tables_to_empty :
	req1 = "SHOW TABLES LIKE " + "'%" + pattern + "%'" + " ;"
	print(req1)
	cur.execute(req1)
	for elt in cur : 
		tables_to_empty2.append(elt)



for elt in tables_to_empty2 :
   req = "TRUNCATE TABLE " + elt[0] +" ;"
   print(req)
   cur.execute(req)

cur.close()

conn.close()

