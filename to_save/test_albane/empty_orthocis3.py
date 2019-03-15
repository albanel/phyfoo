import pymysql

conn = pymysql.connect(host='192.168.1.198',  user='orthocis3', passwd='r77ujtt9', db='orthocis3')
cur = conn.cursor()
req1 = "SHOW TABLES;"
cur.execute(req1)

table_list = []

for table in cur :
   table_list.append(table)

#print(table_list)

for elt in table_list :
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

