from bs4 import BeautifulSoup

with open("database_file.txt","r") as db_html :
	html_doc = db_html.read()
	soup = BeautifulSoup(html_doc)
	p = soup.find_all(id="div_method_link_species_set")
	for elt in p :
		list1 = elt.find_all("td")
	ct = 0
	print(list1[0])
	for elt in list1 :
		if ct % 5 == 0 :
			elt = str(elt)
			elt = elt.replace("<b>","</b>")
			elt = elt.split("</b>")
			print(elt[1])
		ct+=1
	
	
