#!/usr/bin/python
from random import randint
from os import system
nDemes = 200
nInd = 500
nGenerations = 10000
nNtrlLoci = 20
muNtrl = 0.00001
nQuantiLoci = [1, 2, 10, 100]
muQuanti = 0.00001
maxOffs = [3, 10, 100]
immigration = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 2]
extinction = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2]
recolonizer = [1, 5, 10]

for i in nQuantiLoci:
	for j in maxOffs:
		for k in recolonizer:
			for l in immigration:
				for m in extinction:
						seed = randint(1, 1000000)
						commande = "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}".format(nDemes, nInd, \
nGenerations, nNtrlLoci, muNtrl, i, muQuanti, j, l, m, k, seed)
						commande = 'bsub -q normal <<< "quantiSex {0}  | statsQuantiSex.py {0} >>res"'.format(commande)
						print(commande)
						#system(commande)
			#system("sleep 10")
