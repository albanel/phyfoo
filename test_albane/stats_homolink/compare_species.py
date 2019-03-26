non_global_file=open("species_name_to_link.txt","r")
global_file=open("species_name_to_link_global.txt","r")

non_global_dic={}
global_dic={}

whole_dico={}

for line in non_global_file :
	line=line.replace("\n","")
	print(line)
	line=line.split(" ")
	non_global_dic[line[0]]=line[1]

non_global_file.close()

for line2 in global_file :
	print(line2)
	line2=line2.replace("\n","")
	line2=line2.split(" ")
	global_dic[line2[0]]=line2[1]

global_file.close()

print(non_global_dic)

print("=============================================================")

print(global_dic)


for key in non_global_dic.keys() :
	whole_dico[key]=[non_global_dic[key]]
	
for key in global_dic.keys() :
	for key2 in whole_dico.keys() :
		if key==key2 :
			whole_dico[key].append(global_dic[key])


print("=============================================================")

print(whole_dico)
