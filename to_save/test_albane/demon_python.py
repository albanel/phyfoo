#ce script recherche d'abord dans tous les noms de tables puis dans toutes les colonnes le mot clé indiqué en argument
import pymysql
import time
import os

conn = pymysql.connect(host='192.168.1.198',  user='orthocis3', passwd='r77ujtt9', db='orthocis3')
cur = conn.cursor()
req1 = "SELECT COUNT(*) FROM current_global_homolink ;"
#print(req1)

cur.execute(req1)

for response in cur :
        nbre1=response[0]

time.sleep(60)

cur.execute(req1)

for response in cur :
        nbre2=response[0]

while nbre1 != nbre2 :
	nbre1=nbre2
	time.sleep(600)
	cur.execute(req1)
	for response in cur :
		nbre2=response[0]
	print("nbre2",nbre2)

cur.close()
conn.close()

print('\a')
time.sleep(0.5)
print('\a')
time.sleep(0.5)
print('\a')

os.system('echo "test" | mail -s test albane.lysiak@etudiant.univ-rennes1.fr')

