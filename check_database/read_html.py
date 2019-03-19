from bs4 import BeautifulSoup

original_generic_tables = ["genome_db","parameter","member","method_link","method_link_species_set","species_set","homology"]

generic_dico_file = open("generic_dico.txt","w")

genericTableToField = {}

with open("generic_file.txt","r") as db_html :
	html_doc = db_html.read()
	soup = BeautifulSoup(html_doc, features="html.parser")
	for table in original_generic_tables :
		p = soup.find_all(id="div_"+table)
		for elt in p :
			list1 = elt.find_all("td")
		ct = 0
		fields = []
		for elt in list1 :
			if ct % 5 == 0 :
				elt = str(elt)
				elt = elt.replace("<b>","</b>")
				elt = elt.split("</b>")
				fields.append(elt[1])
			ct+=1
		genericTableToField[table] = fields

generic_dico_file.write(str(genericTableToField ))
generic_dico_file.close()
