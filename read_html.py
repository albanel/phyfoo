from bs4 import BeautifulSoup

liste = []

with open("database_file.txt","r") as db_html :
	html_doc = db_html.read()
	soup = BeautifulSoup(html_doc)
	p = soup.find_all(id="div_homology_member")

#print(liste)
