#!/software/bin/python
# -*-coding:Latin-1 -*
#	! /usr/bin/python

import sys
import os
from numpy import array
from numpy import nan
from numpy import nanmean
from numpy import nanstd
from numpy import median

# file = sys.argv[1]
npop = int(sys.argv[1])
nlocus = int(sys.argv[4])

def puissanceDeux(x):   #fonction qui retourne "x" a la puissance deux
        return( x * x )

def expectedHetero(compoAlleliq):       #fonction qui retourne l heterozygotie theorique a partir d'une contingence allelique 
        ntot = float(sum(compoAlleliq.values()))
	if(len(compoAlleliq) == 1):	# Si locus monomorphique, alors Hs = 0
		return(0)
	if(len(compoAlleliq) > 1):	# Si au moins 2 alleles, alors Hs = 1 - (somme des frequences au carre)
       		res = 1
	       	for i in compoAlleliq.values():
			res -= puissanceDeux(( i / ntot ))
	res = res * sum(compoAlleliq.values()) / (sum(compoAlleliq.values()) - 1)
        return(res)

# input = open(file, "r")

nind = []
nTot = 0.
nindTmp = 0.
cnt = 0.
compoPopTot = {}
heterozygotes = {}
nHeteroDeme = {}		#initialiser le dictionnaire contenant le nombre d hetero par locus au sein du nouveau deme
compoDeme = {}		#initialiser le dictionnaire contenant la composition allelique au sein du nouveau deme
HO = {}
HS = {}
femaleAllocation = []

for i in range(nlocus):
	compoPopTot[i] = {}	# un dictionnaire par locus contenant la composition allélique de toute la métapop
	heterozygotes[i] = 0.	# une case par locus contenant le nombre d'hétérozygotes
	nHeteroDeme[i] = 0.
	compoDeme[i] = {}
	HO[i] = []
	HS[i] = []
	
#for i in input:
for i in sys.stdin:
	i = i.strip().split(" ")	#formate la ligne
	population = int(i[1])	#recupere indice de population
	genotypes = i[5:(5+nlocus*2)]	#recupere genotypes
	femaleAllocation.append(float(i[-3]))
	nTot += 1			#compte nombre d'individus dans la metapop
	if cnt!=population:	#si la population actuelle est différente de la precedente alors:
		cnt=population
		nind.append(nindTmp)	#ajouter le nombre d indiv du deme precedent au vecteur de contingence 'nind'
		if nindTmp == 1:
			for k in range(nlocus):
				HS[k].append(nan)
				HO[k].append(nan)
		else:
			for k in range(nlocus):
				HS[k].append(expectedHetero(compoDeme[k]))
				HO[k].append(nHeteroDeme[k]/nindTmp)
		nindTmp = 0.		#initialiser le nombre d'indiv dans les demes 'nindTmp'
		nHeteroDeme = {}		#initialiser le dictionnaire contenant le nombre d hetero par locus au sein du nouveau deme
		compoDeme = {}		#initialiser le dictionnaire contenant la composition allelique au sein du nouveau deme
		for j in range(nlocus):	#réinitialiser les dictionnaires
			nHeteroDeme[j] = 0.
			compoDeme[j] = {}
	nindTmp += 1
	locus = 0.
	for j in range(nlocus):
		alleleA = genotypes[ j * 2 ]
		alleleB = genotypes[ j * 2 + 1 ]
		if alleleA != alleleB:
			heterozygotes[locus] += 1
			nHeteroDeme[locus] += 1
		if alleleA in compoPopTot[locus] and alleleA in compoDeme[locus] :
			compoPopTot[locus][alleleA] += 1
			compoDeme[locus][alleleA] += 1
		if alleleA in compoPopTot[locus] and alleleA not in compoDeme[locus] :
			compoPopTot[locus][alleleA] += 1
			compoDeme[locus][alleleA] = 1
		if alleleA not in compoPopTot[locus] :
			compoPopTot[locus][alleleA] = 1
			compoDeme[locus][alleleA] = 1
		if alleleB in compoPopTot[locus] and alleleB in compoDeme[locus] :
			compoPopTot[locus][alleleB] += 1
			compoDeme[locus][alleleB] += 1
		if alleleB in compoPopTot[locus] and alleleB not in compoDeme[locus] :
			compoPopTot[locus][alleleB] += 1
			compoDeme[locus][alleleB] = 1
		if alleleB not in compoPopTot[locus] :
			compoPopTot[locus][alleleB] = 1
			compoDeme[locus][alleleB] = 1
		locus += 1
#input.close()

#ajoute stats du dernier deme
nind.append(nindTmp)
if nindTmp == 1:
	for k in range(nlocus):
		HS[k].append(nan)
		HO[k].append(nan)
else:
	for k in range(nlocus):
		HS[k].append(expectedHetero(compoDeme[k]))
		HO[k].append(nHeteroDeme[k]/nindTmp)

HS_tot = {}
HO_tot = {}
HT_tot = {}
HS_tot2 = []
HO_tot2 = []
HT_tot2 = []
Fst = []
Fis = []
Fit = []
DJost = []

for i in range(nlocus):
	HS_tmp = array(HS[i])
	HO_tmp = array(HO[i])
	HS_tot[i] = nanmean(HS_tmp)	
	HO_tot[i] = nanmean(HO_tmp)
	HT_tot[i] = expectedHetero(compoPopTot[i]) 
	HS_tot2.append(HS[i])
	HO_tot2.append(HO[i])
	HT_tot2.append(expectedHetero(compoPopTot[i]))
#	fo rj in range(len(nind)):	# boucle sur les demes (nind est un vecteur contenant le #_of_individuals par deme)
#		if nind[j]>1 and HS[i][j]>0 and HO[i][j]>0:
#			HS_tot[i] += HS[i][j]#*nind[j]
#			HO_tot[i] += HO[i][j]#*nind[j]
#			tmp += nind[j]
#			tmp += 1
#	if tmp == 0:
#		HS_tot[i] = "NA"
#		HO_tot[i] = "NA"
#		HT_tot[i] = expectedHetero(compoPopTot[i])
#	if tmp > 0:
#		HS_tot[i] /= tmp
#		HO_tot[i] /= tmp
#		HT_tot[i] = expectedHetero(compoPopTot[i])
HS_tot2 = array(HS_tot2)
HO_tot2 = array(HO_tot2)
HT_tot2 = array(HT_tot2)

# calcul du Fst
for i in range(nlocus):
	if HS_tot[i] == nan or HT_tot[i] == 0 or HT_tot[i] == nan:
		Fst.append(0)
	else:
		Fst.append(round(1-HS_tot[i] / HT_tot[i], 5))

# calcul du Fis
for i in range(nlocus):
	if HO_tot[i] == nan or HS_tot[i] == 0 or HS_tot[i] == nan:
		Fis.append(0)
	else:
		Fis.append(round(1-HO_tot[i] / HS_tot[i], 5))

# calcul du Fit
for i in range(nlocus):
	if HO_tot[i] == nan or HT_tot[i] == 0 or HT_tot[i] == nan:
		Fit.append(0)
	else:
		Fit.append(round(1-HO_tot[i]/ HT_tot[i], 5))

# calcul du DJost
for i in range(nlocus):
	if HT_tot[i] == nan or HS_tot[i] == nan:
		DJost.append(nan)
	else:
		DJost.append( round(((HT_tot[i]-HS_tot[i]) / (1.-HS_tot[i])) * (npop) / (npop-1.) , 5))

# mettre des "NA" pour les locus avec des Fst négatifs
for i in range(len(Fst)):
	if(Fst[i] < 0):
		Fst[i] = nan
		Fis[i] = nan
		Fit[i] = nan
		DJost[i] = nan

Fst = array(Fst)
Fis = array(Fis)
Fit = array(Fit)
DJost = array(DJost)
femaleAllocation = array(femaleAllocation)

#output="I\tE\tnTot\tAllocF\tAllocF_std\tFst_CR\tFst_SN\tFst_WC\tFst_std\tFis_CR\tFis_SN\tFis_WC\tFis_std\tFit_CR\tFit_SN\tFit_WC\tFit_std\tDJost\n{0}\t{1}\t{17}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\n".format(I, E, allocation, allocation_std, median(Fst), FstQuanti, FstQuantiWC, FstQuantiWC_std, median(Fis), FisQuanti, FisQuantiWC, FisQuantiWC_std, median(Fit), FitQuanti, FitQuantiWC, FitQuantiWC_std, median(DJost), nTot)
output1 = "nDemes\tnIndPerDem\tnGenerations\tnNtrlLoci\tmuNtrl\tnQuantiLoci\tmuQuanti\tmaxNoffspring\timmigration\textinction\tnRecolonizer\tseed\tFst_median\tFst_mean\tFst_std\tFis_median\tFis_mean\tFis_std\tFit_median\tFit_mean\tFit_std\tDJost_median\tDJost_mean\tDJost_std\tFemaleAlloc_median\tFemaleAlloc_mean\tFemaleAlloc_std\tHS_median\tHS_mean\tHS_std\tHO_median\tHO_mean\tHO_std\tHT_median\tHT_mean\tHT_std\n"
#for i in range(nlocus):
#	output += "{0}\t{1}\t{2}\t{3}\n".format(Fst[i], Fis[i], Fit[i], DJost[i])
output2 = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{25}\t{26}\t{27}\t{28}\t{29}\t{30}\t{31}\t{32}\t{33}\t{34}\t{35}".format(\
sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],  sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], \
round(median(Fst), 5), round(nanmean(Fst), 5), round(nanstd(Fst), 5), \
round(median(Fis), 5), round(nanmean(Fis), 5), round(nanstd(Fis), 5), \
round(median(Fit), 5), round(nanmean(Fit), 5), round(nanstd(Fit), 5), \
round(median(DJost), 5), round(nanmean(DJost), 5), round(nanstd(DJost), 5), \
round(median(femaleAllocation), 5), round(nanmean(femaleAllocation), 5), round(nanstd(femaleAllocation), 5), \
round(median(HS_tot2), 5), round(nanmean(HS_tot2), 5), round(nanstd(HS_tot2), 5), \
round(median(HO_tot2), 5), round(nanmean(HO_tot2), 5), round(nanstd(HO_tot2), 5), \
round(median(HT_tot2), 5), round(nanmean(HT_tot2), 5), round(nanstd(HT_tot2), 5), \
)

print(output2)
outputFile=open(sys.argv[1] + sys.argv[2] + sys.argv[3] + sys.argv[4] + sys.argv[5] + sys.argv[6] + sys.argv[7] + sys.argv[8] + sys.argv[9] + sys.argv[10] + sys.argv[11] + sys.argv[12] + ".txt", "w")
outputFile.write(output1+output2)
outputFile.close()

