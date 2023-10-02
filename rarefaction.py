import sys
import matplotlib.pyplot as plt
import random
from statistics import median

orthofinder_dir = sys.argv[0]
if orthofinder_dir[-1] == "/":
	orthofinder_dir = orthofinder_dir[:-1]
orthofile = "Orthogroups.GeneCount.tsv"

figname = sys.argv[-1]
spplist = []

with open(orthofile) as file:
	lines = file.readlines()
	spplist = lines[0].strip().split("\t")[1:-1]

spppos = {i: name for i, name in enumerate(spplist)}

setdict = {name: set() for name in spplist}

for line in lines[1:]:
	data = line.strip().split("\t")
	orthogroup = data[0]
	values = [int(x) if x.isdigit() else 0 for x in data[1:-1]]
	# values = []
	# for x in data[1:-1]:
	# 	if x.isdigit():
	# 		values.append(int(x))
	# 	else:
	# 		print(f"Error in Orthogroup: {orthogroup} ,  {x} is not a digit")
	# 		values.append(0)
	for i, value in enumerate(values):
		if value > 0:
			setdict[spplist[i]].add(orthogroup)

rand_pan, rand_core = [], []

for i in spplist:
	rand_pan.append([])
	rand_core.append([])

while len(rand_pan[0]) < 1000:
	names, pan, core = [], [], []
	Aset, Bset = set(), set()
	random.shuffle(spplist)

	for spp in spplist:
		names.append(spp)
		if len(pan) == 0:
			pan.append(len(setdict[spp]))
			core.append(len(setdict[spp]))
			Aset, Bset = setdict[spp], setdict[spp]
		else:
			Aset, Bset = Aset | setdict[spp], Bset & setdict[spp]
			pan.append(len(Aset))
			core.append(len(Bset))

	for i in range(0, len(pan)):
		rand_pan[i].append(pan[i])
		rand_core[i].append(core[i])

med_pan, med_core = [], []

for i in rand_pan:
	med_pan.append(median(i))
for i in rand_core:
	med_core.append(median(i))

fig = plt.figure()
p1 = plt.violinplot(rand_pan[1:], positions=range(2,len(spplist)+1))
p2 = plt.violinplot(rand_core)
p3 = plt.plot(range(1,len(spplist)+1), med_pan, "b-")
p4 = plt.plot(range(1,len(spplist)+1), med_core, "r-")
plt.axhline(y=med_pan[-1], color='#CCCCCC', linestyle='-.')
plt.axhline(y=med_core[-1], color='#CCCCCC', linestyle='-.')
plt.xrange=[1,len(spplist)]
plt.text(len(spplist), med_pan[-1]+20, str(int(med_pan[-1])), fontsize=12)
plt.text(len(spplist), med_core[-1]+20, str(int(med_core[-1])), fontsize=12)
plt.title(figname)

#p2 = plt.bar(x=range(1,len(names)+1), height=core, color="tab:blue")
#plt.xticks(rotation = 90)
plt.ylabel('Cumulative number of orthogroups', fontsize=12)

fig.set_size_inches(18.5, 10.5)
plt.savefig("rarefactioncurve_"+figname+".svg")
plt.savefig("rarefactioncurve_"+figname+".png")
plt.show()