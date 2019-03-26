import matplotlib.pyplot as plt

x=[]
#y=[]

with open("global_homolink_distribution.txt","r") as homolink :
	for line in homolink :
		line=line.split("\t")
		for nbre in line :
			if nbre != "" : 
				x.append(int(nbre))
		#y.append(int(nbre_liens))		


#print(x)

#print(y)

# plotting a histogram 

plt.hist(x, bins=30)
  
# naming the x axis 
plt.ylabel("Nombre de gènes") 
# naming the y axis 
plt.xlabel("Connectivité") 
  
# giving a title to my graph 
plt.title("Connectivité des gènes dans la table global_homolink") 
  
# save the plot
plt.savefig("global_stats.png")
