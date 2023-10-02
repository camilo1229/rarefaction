import sys
import matplotlib.pyplot as plt
import random
from statistics import median, mean

# Rutas de los archivos TSV de OrthoFinder y listas de especies
orthofinder_dir = sys.argv[0]
if orthofinder_dir[-1] == "/":
    orthofinder_dir = orthofinder_dir[:-1]

orthofile = "Orthogroups.GeneCount.tsv"  # Reemplaza con la ubicación real de tu archivo TSV
sppfile1 = "Orthogroups.GeneCount.tsv"  # Reemplaza con la ubicación real del archivo de cultivos
sppfile2 = "Orthogroups.GeneCount.tsv"  # Reemplaza con la ubicación real del archivo silvestres

figname = sys.argv[-1]
spppos = {}
firstlist, secondlist = [], []

# Leer la lista de especies de cultivos
with open(sppfile1) as file:
    firstlist = [line.strip() for line in file]

# Leer la lista de especies silvestres
with open(sppfile2) as file:
    secondlist = [line.strip() for line in file]

# Combinar ambas listas de especies
spplist = firstlist + secondlist

# Crear un diccionario para almacenar los conjuntos de Orthogroups por especie
setdict = {name: set() for name in spplist}

# Leer el archivo TSV de OrthoFinder
with open(orthofile) as file:
    for line in file:
        chunk = line.strip().split("\t")
        if chunk[0] == "Orthogroup":
            for i in range(len(chunk)):
                if chunk[i] in spplist:
                    spppos[i] = chunk[i]
                    setdict[chunk[i]] = set()
        else:
            for i in range(len(chunk)):
                if i in spppos and i != 0 and i != len(chunk) - 1 and int(chunk[i]) > 0:
                    setdict[spppos[i]].add(chunk[0])

rand_pan, rand_core = [], []

for i in spplist:
    rand_pan.append([])
    rand_core.append([])

while len(rand_pan[0]) < 1000:
    names, pan, core = [], [], []
    Aset, Bset = set(), set()
    random.shuffle(firstlist)
    random.shuffle(secondlist)
    newlist = firstlist + secondlist
    for spp in newlist:
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
    med_pan.append(mean(i))
for i in rand_core:
    med_core.append(mean(i))

fig = plt.figure()
plt.axhline(y=med_pan[-1], color='#CCCCCC', linestyle='-.')
plt.axhline(y=med_core[-1], color='#CCCCCC', linestyle='-.')
plt.axvline(x=len(firstlist), color='#CCCCCC', linestyle='-.')
p1 = plt.violinplot(rand_pan[1:], positions=range(2, len(spplist) + 1))
p2 = plt.violinplot(rand_core)
p3 = plt.plot(range(1, len(spplist) + 1), med_pan, "b-")
p4 = plt.plot(range(1, len(spplist) + 1), med_core, "r-")

plt.text(len(firstlist) - 5, med_pan[-1] + 20, str("Cultivars"), fontsize=12)
plt.text(len(firstlist) + 2, med_pan[-1] + 20, str("Wild"), fontsize=12)

plt.xlim(1, len(spplist))
plt.text(len(spplist), med_pan[-1] + 20, str(int(med_pan[-1])), fontsize=12)
plt.text(len(spplist), med_core[-1] + 20, str(int(med_core[-1])), fontsize=12)
plt.text(len(firstlist) - 3, med_pan[len(firstlist) - 1] + 100, str(int(med_pan[len(firstlist) - 1])), fontsize=12)
plt.text(len(firstlist), med_core[len(firstlist) - 1] + 100, str(int(med_core[len(firstlist) - 1])), fontsize=12)
plt.title(figname)

plt.ylabel('Número acumulativo de grupos ortológicos', fontsize=12)

fig.set_size_inches(18.5, 10.5)
plt.savefig("rarefactioncurve_" + figname + ".svg")
plt.savefig("rarefactioncurve_" + figname + ".png")
plt.show()
