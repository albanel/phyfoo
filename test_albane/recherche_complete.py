#ce script recherche d'abord dans tous les noms de tables puis dans toutes les colonnes le mot clé indiqué en argument
import pymysql
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("mot")
args = parser.parse_args()
conn = pymysql.connect(host='192.168.1.198',  user='orthocis3', passwd='r77ujtt9', db='orthocis3')
cur = conn.cursor()
req1 = "SHOW TABLES FROM orthocis3 LIKE '%"+args.mot+"%' ;"
#print(req1)
cur.execute(req1)
for response in cur:
   print("table",str(response))

req2 = "SELECT column_name, table_name FROM information_schema.columns WHERE column_name LIKE '%"+args.mot+"%' ;"

#print(req2)
cur.execute(req2)
for response in cur:
   print("column",str(response))

cur.close()

conn.close()

