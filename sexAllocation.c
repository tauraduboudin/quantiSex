#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_permutation.h>
#define VERSION "09.05.2016"
#define DEPENDENCY "diveRsity.R\n"
#define MAX_NUMBER_OF_INITIAL_NTRL_ALLELES 999	// number of segregating alleles when generating the first parental population
#define RANGE 0.1	// value in [0;1] to modify the current allelic effect between [(1-RANGE) x current_value ; (1+RANGE) * current_value].
#define KRED  "\033[1m\033[31m"
#define KMAG  "\x1B[31m"
#define STOP  "\x1B[0m"

//	gcc sexAllocation.c -L/usr/local/lib -lgsl -lgslcblas -lm -Wall -Wextra -Wshadow -Werror -O3 -o quantiSex
//	./quantiSex 100 200 100 10 0.00001 1 0.0001 10 0 0.3 1 0 1 123

typedef struct Deme Deme;
struct Deme{
	int nIndividus;
	long* ntrlLoci;	// contains the neutral alleles for nIndividus (number of individuals within the deme) x 2 (diploids) x nNtlrLoci (number of neutral loci)
	double* quantiLoci; // contains the allelic effects for nIndividus (number of individuals within the deme) x 2 (diploids) x nQuantiLoci (number of quantitative loci)
	int* sexChro; // two values per individual: 00 = XX (or ZZ), 01 = XY (or ZW)
	int* sex; // one value per individual: 0 = heterogametic, 1 = homogametic 
	double* femaleAllocation; // =sum of the allelic effects over the quantitative loci
	double* maleAllocation; // =(1 - femaleAllocation)
	int* nOffsprings; // floor(X) + Binom(n=1, p=X-floor(x)) where X = femaleAllocation x fecundity
};

void initializePopulation(gsl_rng *r, Deme* population, const int nDemes, const int maxIndPerDem, const int nNtrlLoci, const int nQuantiLoci, const int fecundity, const double sexAvantage, const int sexualSystem);
void libererMemoirePopulation(Deme* population, const int nDemes);
void afficherPopulation(Deme* population, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const int sexualSystem);
void configMetapop(gsl_rng* r, Deme* population, const int nDemes, const double migration, const double extinction, int nImmigrants[], int extinctionStatus[], int nProducedSeeds[], const int maxIndPerDem);
void setToZero(const int nDemes, int nImmigrants[], int extinctionStatus[], int nProducedSeeds[]);
void initializeNewPopulation(Deme* newPopulation, const int nDemes, const int maxIndPerDem, const int nProducedSeeds[], const int nNtrlLoci, const int nQuantiLoci, const int fecundity);
void panmixie(gsl_rng* r, Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const double ntrlMutation, const double quantiMutation, const int fecundity, const double sexAvantage, const int sexualSystem);
void weightedSample(gsl_rng* r, const double* liste, const double* weights, double* target, const int sizeOfListe, const int nTrials);
void replacement(gsl_rng* r, Deme* population, Deme* newPopulation, const int nDemes, const int maxIndPerDem, const int nNtrlLoci, const int nQuantiLoci, const int nImmigrants[], int nProducedSeeds[], int extinctionStatus[], const int recolonization, const int generation, const int sexualSystem);
void writeNindividuals(const Deme* population, const int nDemes, const double extinction, const double migration, const int seed);
void genePop(Deme* population, const int nDemes, const int nNtrlLoci, const int seed, int time);
void checkCommandLine(int argc);
void statisticsPopulations(Deme* population, const int nDemes, const int maxIndPerDem, const int nQuantiLoci, const int fecundity, const double migration, const double extinction, const int recolonization, const int sexualSystem, const double sexAvantage, const int seed, int time);
double fst(const int maxIndPerDem, const double extinction, const int recolonization, const double migration);

int main(int argc, char *argv[]){

	checkCommandLine(argc); // stop the code if the number of arguments doesn't fit with the expected one

	int i = 0;
	// Get Parameters from comamnd line
	const int nDemes = atoi(argv[1]); // number of demes
	const int maxIndPerDem = atoi(argv[2]);	// carrying capacity per deme
	const int nGeneration = atoi(argv[3]);	// number of generations to simulate
	const int nNtrlLoci = atoi(argv[4]);	// number of neutral loci
	const double ntrlMutation = atof(argv[5]);	// mutation rate of the ntrl loci
	const int nQuantiLoci = atoi(argv[6]);	// number of quantitative loci
	const double quantiMutation = atof(argv[7]);	// mutation rate of the quantative loci
	const int fecundity = atoi(argv[8]);	// max number of offspring when femAlloc=100%
	const double migration = atof(argv[9]);	// immigration rate
	const double extinction = atof(argv[10]);	// extinction rate
	const int recolonization = atoi(argv[11]);	// number of recolonizing individuals
	const int sexualSystem = atoi(argv[12]);    // 0 = only hermaphrodites; 1 = XY system; 2 = ZW system
	const double sexAvantage = atof(argv[13]); // avantage confered by the Y or Z chromosome over hermaphrodites
	const int seed = atoi(argv[14]);

	// Random generator
        const gsl_rng_type *T;
        gsl_rng *r;
        gsl_rng_env_setup();

        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, seed);

	// Initializing the metapopulation
	Deme* population = NULL;
	population = malloc(nDemes * sizeof(Deme));
	if(population == NULL){
		exit(0);
	}

	initializePopulation(r, population, nDemes, maxIndPerDem, nNtrlLoci, nQuantiLoci, fecundity, sexAvantage, sexualSystem);

	// Evolution of the metapopulation
	for(i=0; i<=nGeneration; i++){	// start of the loop 'i' over the generations

		int* nImmigrants = NULL; // stock the number of received immigrants per deme
		int* extinctionStatus = NULL;	// per deme: 0 = non-extincted; 1 = extincted
		int* nProducedSeeds = NULL; // stock the numer of produced seeds per deme
		Deme* newPopulation = NULL;
		nImmigrants = malloc(nDemes * sizeof(int));
		extinctionStatus = malloc(nDemes * sizeof(int));
		nProducedSeeds = malloc(nDemes * sizeof(int));
		newPopulation = malloc(nDemes * sizeof(Deme));

		if(nImmigrants == NULL || extinctionStatus == NULL || nProducedSeeds == NULL || newPopulation == NULL){
			exit(0);
		}
		setToZero(nDemes, nImmigrants, extinctionStatus, nProducedSeeds);	// initialize vectors used to configure the new-population

		configMetapop(r, population, nDemes, migration, extinction, nImmigrants, extinctionStatus, nProducedSeeds, maxIndPerDem);	// get the parameters to create the new-population

		initializeNewPopulation(newPopulation, nDemes, maxIndPerDem, nProducedSeeds, nNtrlLoci, nQuantiLoci, fecundity);

		panmixie(r, population, newPopulation, nDemes, nNtrlLoci, nQuantiLoci,  ntrlMutation, quantiMutation, fecundity, sexAvantage, sexualSystem);

		statisticsPopulations(newPopulation, nDemes, maxIndPerDem, nQuantiLoci, fecundity, migration, extinction, recolonization, sexualSystem, sexAvantage, seed, i);
		if(i == nGeneration){
			genePop(newPopulation, nDemes, nNtrlLoci, seed, i);
		}
		//writeNindividuals(newPopulation, nDemes, extinction, migration, seed); // to un-comment only if we want #individuals and femaleAllocation per deme per generation

		// remplacer population par newPopulation
		// Nettoyer memoire
		libererMemoirePopulation(population, nDemes);
		free(population);
		population = malloc(nDemes * sizeof(Deme));

		if(population == NULL){
			exit(0);
		}

		replacement(r, population, newPopulation, nDemes, maxIndPerDem, nNtrlLoci, nQuantiLoci, nImmigrants, nProducedSeeds, extinctionStatus, recolonization, i, sexualSystem); // replace the parents (poplation) by the offspring (newPopulation)

		free(nImmigrants);
		free(extinctionStatus);
		free(nProducedSeeds);
		libererMemoirePopulation(newPopulation, nDemes);
		free(newPopulation);
	}	// end of the loop 'i' over the generations
	libererMemoirePopulation(population, nDemes);
	free(population);
	return(0);
}


void initializePopulation(gsl_rng* r, Deme* population, const int nDemes, const int maxIndPerDem, const int nNtrlLoci, const int nQuantiLoci, const int fecundity, const double sexAvantage, const int sexualSystem){
	int i = 0;
	int j = 0;
	int k = 0;
	int valuesNtrlAlleles = MAX_NUMBER_OF_INITIAL_NTRL_ALLELES;
	double minQuanti = 0.0;
	double maxQuanti = 0.0;
	maxQuanti = 1.0/2/nQuantiLoci;
//	maxQuanti = 0.5/2/nQuantiLoci;	// uncomment to fix sex allocation to 0.5
	for(i=0; i<nDemes; i++){	// loop along the demes
		population[i].nIndividus = maxIndPerDem;
		population[i].ntrlLoci = malloc(2 * maxIndPerDem * nNtrlLoci * sizeof(long));
		population[i].quantiLoci = malloc(2 * maxIndPerDem * nQuantiLoci * sizeof(long));
		population[i].sexChro = malloc(2 * maxIndPerDem * sizeof(int));
		population[i].sex = malloc(maxIndPerDem * sizeof(int));
		population[i].femaleAllocation = malloc(maxIndPerDem * sizeof(long));
		population[i].maleAllocation = malloc(maxIndPerDem * sizeof(long));
		population[i].nOffsprings = malloc(maxIndPerDem * sizeof(int));

		if(population[i].ntrlLoci == NULL || population[i].quantiLoci == NULL || population[i].sexChro == NULL || population[i].sex == NULL || population[i].femaleAllocation == NULL || population[i].maleAllocation == NULL || population[i].nOffsprings == NULL){
			exit(0);
		}

		int cnt = -1;
		for(j=0; j<maxIndPerDem; j++){	// loop along the individuals
			cnt += 1;
			population[i].femaleAllocation[j] = 0.0;
			population[i].maleAllocation[j] = 1.0;
			for(k=0; k<(2*nNtrlLoci); k++){	// loop along the {2: diploid} x {nNtrlLoci: number of neutral loci} positions
				population[i].ntrlLoci[j*2*nNtrlLoci+k] = gsl_rng_uniform_int(r, valuesNtrlAlleles) + 1;
			}
			for(k=0; k<(2*nQuantiLoci); k++){	// loop along the {2: diploid} x {nQuantiLoci: number of quantitative loci} positions
				population[i].quantiLoci[j*2*nQuantiLoci+k] = gsl_ran_flat(r, minQuanti, maxQuanti);
//				population[i].quantiLoci[j*2*nQuantiLoci+k] = maxQuanti; // uncomment to fix sex allocation to 0.5
				population[i].femaleAllocation[j] += population[i].quantiLoci[j*2*nQuantiLoci+k];
				population[i].maleAllocation[j] -= population[i].quantiLoci[j*2*nQuantiLoci+k];
			}
			if(sexualSystem == 0){
				population[i].sexChro[2*j] = 0;
				population[i].sexChro[2*j + 1] = 0;
				population[i].sex[j] = 1;
			}else{
				if(cnt%2 == 0){	// homogametic individual
					population[i].sexChro[2*j] = 0;
					population[i].sexChro[2*j + 1] = 0;
					population[i].sex[j] = 1;
				}
				if(cnt%2 != 0){ // heterogametic individual
					population[i].sexChro[2*j] = 0;
					population[i].sexChro[2*j + 1] = 1;
					population[i].sex[j] = 0;
				}
			}

			//printf("%d\n%d\n", population[i].sexChro[2*j], population[i].sexChro[2*j+1]);
            population[i].nOffsprings[j] = floor(fecundity * population[i].femaleAllocation[j]) + gsl_ran_binomial(r, (fecundity * population[i].femaleAllocation[j]) - floor(fecundity * population[i].femaleAllocation[j]), 1);	// nOffs = floor(fecundity x femaleAllocation) + 1 according to a random Binomial integer

            if(population[i].sexChro[2*j] != population[i].sexChro[2*j + 1]){ // if heterogametic individual
                if(sexualSystem == 1){  // if XY system
                    population[i].femaleAllocation[j] = 0;
                    population[i].maleAllocation[j] = sexAvantage;
                    population[i].nOffsprings[j] = 0;
                }
                if(sexualSystem == 2){  // if ZW system
                    population[i].femaleAllocation[j] = sexAvantage;
                    population[i].maleAllocation[j] = 0;
                    population[i].nOffsprings[j] = floor(fecundity * sexAvantage) + gsl_ran_binomial(r, (fecundity * sexAvantage) - floor(fecundity * sexAvantage), 1);
                }
            }
        }	// end of loop along the individuals
    }	// end of loop along the demes
}

void setToZero(const int nDemes, int nImmigrants[], int extinctionStatus[], int nProducedSeeds[]){
	// function to set all to demes to zero for the #of immigrants, the extinction status and the #of produced seeds.
	// before configMetapop()
	int i = 0;
	for(i=0; i<nDemes; i++){
		nImmigrants[i] = 0;
		extinctionStatus[i] = 0;
		nProducedSeeds[i] = 0;
	}
}

void configMetapop(gsl_rng* r, Deme* population, const int nDemes, const double migration, const double extinction, int nImmigrants[], int extinctionStatus[], int nProducedSeeds[], const int maxIndPerDem){
	// function to get per deme the #of immigrants received, the extinction status and the #of produced seeds.
	// after setToZero()
	// modifies nImmigrants[]; extinctionStatus[] and nProducedSeeds[]
	int i = 0;
	int j = 0;
	const unsigned int binomTrials = 1;

	for(i=0; i<nDemes; i++){	// start of the loop over the demes
		int nHeterogametic = 0;  // number of heterogametic individuals in the deme

		for(j=0; j<population[i].nIndividus; j++){
			if(population[i].sex[j] == 0){
				nHeterogametic += 1;
			}
		}
		if(nHeterogametic == population[i].nIndividus){
			printf("PROUT deme: %d\n", i);
			extinctionStatus[i] = 1;
		}

		nImmigrants[i] = gsl_ran_poisson(r, migration);
		if(nImmigrants[i] > maxIndPerDem){
			nImmigrants[i] = maxIndPerDem;
		}
		extinctionStatus[i] = gsl_ran_binomial(r, extinction, binomTrials);	// 0: non-extincted; 1: extincted

		if(extinctionStatus[i] == 1){	// comment this block if you allow recolonization and migration occuring at the same time
			nImmigrants[i] = 0;	// no migrant if a deme is extinct
		}
		for(j=0; j<population[i].nIndividus; j++){	// start of the loop over the individuals
			nProducedSeeds[i] += population[i].nOffsprings[j];
		}	// end of the loop over the individuals
	}	// end of the loop over the demes
}

void initializeNewPopulation(Deme* newPopulation, const int nDemes, const int maxIndPerDem, const int nProducedSeeds[], const int nNtrlLoci, const int nQuantiLoci, const int fecundity){
	// function to initialize the new population: allocate memory and set values to 0
	int i = 0;
	int j = 0;
	int k = 0;
	for(i=0; i<nDemes; i++){ // start of the loop along demes
		int taille = 0;
		taille = nProducedSeeds[i];	// size of the deme = nProducedSeeds within the parental population
		if(taille < fecundity){	// if not enough seeds are produced ==> deme is considered as extincted
			taille = fecundity;
		}
		if(taille > maxIndPerDem){	// if too many individuals have to be present in the deme ==> cutoff to the maxIndPerDem (=carrying capacity)
			taille = maxIndPerDem;
		}

		newPopulation[i].nIndividus = taille;

		newPopulation[i].ntrlLoci = NULL;
		newPopulation[i].quantiLoci = NULL;
		newPopulation[i].sexChro = NULL;
		newPopulation[i].sex = NULL;
		newPopulation[i].femaleAllocation = NULL;
		newPopulation[i].maleAllocation = NULL;
		newPopulation[i].nOffsprings = NULL;

		newPopulation[i].ntrlLoci = malloc(2 * taille * nNtrlLoci * sizeof(long));
		newPopulation[i].quantiLoci = malloc(2 * taille * nQuantiLoci * sizeof(long));
		newPopulation[i].sexChro = malloc(2 * maxIndPerDem * sizeof(int));
		newPopulation[i].sex = malloc(maxIndPerDem * sizeof(int));
		newPopulation[i].femaleAllocation = malloc(taille * sizeof(long));
		newPopulation[i].maleAllocation = malloc(taille * sizeof(long));
		newPopulation[i].nOffsprings = malloc(taille * sizeof(int));

		if(newPopulation[i].ntrlLoci == NULL || newPopulation[i].quantiLoci == NULL || newPopulation[i].sexChro == NULL || newPopulation[i].sex == NULL || newPopulation[i].femaleAllocation == NULL || newPopulation[i].maleAllocation == NULL || newPopulation[i].nOffsprings == NULL){
			exit(0);
		}

		for(j=0; j<taille; j++){ // start the loop along individuals
			newPopulation[i].femaleAllocation[j] = 0.0;
			newPopulation[i].maleAllocation[j] = 0.0;
			newPopulation[i].nOffsprings[j] = 0;

			newPopulation[i].sexChro[2*j] = 0;
			newPopulation[i].sexChro[2*j + 1] = 0;
			newPopulation[i].sex[j] = 0;

			for(k=0; k<(2*nNtrlLoci); k++){ // loop along the {2: diploid} x {nNtrlLoci: number of neutral loci} positions
				newPopulation[i].ntrlLoci[j*2*nNtrlLoci + k] = 0;
			}
			for(k=0; k<(2*nQuantiLoci); k++){       // loop along the {2: diploid} x {nQuantiLoci: number of quantitative loci} positions
				newPopulation[i].quantiLoci[j*2*nQuantiLoci + k] = 0;
			}
		}	// end of the loop along individuals
	} // end of the loop along demes
}

void panmixie(gsl_rng* r, Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const double ntrlMutation, const double quantiMutation, const int fecundity, const double sexAvantage, const int sexualSystem){
	// function returning a new deme after a run of panmixia from the old deme
	// here: only "deme specific" events are simulated (meiosis + mutation)
	double currentAllelicEffect = 0.0;	// current allelic effect of a quantitative allele
	double minEffect = 0.0;	// minimum value that a quantitative mutation can take
	double maxEffect = 0.0;	// maximum value that a quantitative mutation can take
	int i = 0;
	int j = 0;
	int tmp = 0;
	int nMutation = 0;
		for(i=0; i<nDemes; i++){	// start the loop over the nDemes
			int K = 0;	// #of_individuals in the parental deme.
			int N = 0;	// #of_babies to produce
			double* parentalIndexes = NULL; // contain de the parental indexes (i.e:{0, 1, 2, 3, 4} if K == 5)
			double* mothers = NULL;	// array containing the mothers ID. Ex: {2, 19, 3, 3, 1} means that the first baby has individual#2 as mother, second baby has individual#19 as mother, babies 3 and 4 have individual#3
			double* fathers = NULL;	// array containing the fathers ID.

			K = population[i].nIndividus;	// #of_individuals in the parental deme
			N = newPopulation[i].nIndividus; 	// #of_produced_seeds by the K parents

			parentalIndexes = malloc(K * sizeof(double));
			mothers = malloc(N * sizeof(double));	// indexes of mothers of the N autochtones within the deme of size newPopulation[i].nIndividus
			fathers = malloc(N * sizeof(double));

			if(parentalIndexes == NULL || mothers == NULL || fathers == NULL){
				exit(0);
			}

			for(j=0; j<K; j++){
				parentalIndexes[j] = j;
			}
			for(j=0; j<N; j++){
				mothers[j] = 0;
				fathers[j] = 0;
			}

			// Sampling the parents
			weightedSample(r, parentalIndexes, population[i].femaleAllocation, mothers, K, N); // return the indexes of N mothers (in 'mothers') that will be use to generate the N babies
			weightedSample(r, parentalIndexes, population[i].maleAllocation, fathers, K, N); // return the indexes of N fathers (in 'fathers') that will be use to generate the N babies

			// Meiosis and transmission of gametes
			int pos = 0;	// position in deme.ntrlLoci (or deme.quantiLoci) of offsprings
			int pos2 = 0;	// position in deme.ntrlLoci (or deme.quantiLoci) of parents
			for(j=0; j<N; j++){	// loop over the babies.  Transmission of alleles, meiosis = gsl_ran_binomial
				for(tmp=0; tmp<nNtrlLoci; tmp++){	// loop along the neutral loci
					pos = j*nNtrlLoci*2 + tmp*2 + 0;	// position where to put the allele from the mother
					pos2 = mothers[j]*nNtrlLoci*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1); // binomial transmission of the parental alleles with p=0.5 (=recombination at meiosis)
					newPopulation[i].ntrlLoci[pos] = population[i].ntrlLoci[pos2];

					pos += 1;	// position where to put the allele from the father = position of mother_allele + 1
					pos2 = fathers[j]*nNtrlLoci*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1);
					newPopulation[i].ntrlLoci[pos] = population[i].ntrlLoci[pos2];
				}	// end of loop along the neutral loci

				for(tmp=0; tmp<nQuantiLoci; tmp++){	// loop along the quantitative loci
					pos = j*nQuantiLoci*2 + tmp*2 + 0;
					pos2 = mothers[j]*nQuantiLoci*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1);
					newPopulation[i].quantiLoci[pos] = population[i].quantiLoci[pos2];

					pos += 1;
					pos2 = fathers[j]*nQuantiLoci*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1);
					newPopulation[i].quantiLoci[pos] = population[i].quantiLoci[pos2];
				}	// end of loop along the quantitative loci
				
				for(tmp=0; tmp<1; tmp++){ // transmission of the sex chromosomes
					pos = j*2 + tmp*2 + 0;
	               			pos2 = mothers[j]*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1); // sex chromosome from the mother
			                newPopulation[i].sexChro[pos] = population[i].sexChro[pos2];

					pos += 1;
                			pos2 = fathers[j]*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1); // sex chromosome from the father
		        	        newPopulation[i].sexChro[pos] = population[i].sexChro[pos2];
				}

			}	// end of loop along the individuals

			for(j=0; j<N; j++){	// loop over babies. Determine their sex 0 (heterogametic) or 1 (homogametic)
				if(newPopulation[i].sexChro[2*j] == newPopulation[i].sexChro[2*j + 1]){ // if homogametic, sex = 1
					newPopulation[i].sex[j] = 1;
				}else{
					newPopulation[i].sex[j] = 0; // if heterogametic, sex = 0
				}
			}

			// Put mutations
			nMutation = 0;
			nMutation = gsl_ran_binomial(r, ntrlMutation, 2*N*nNtrlLoci);
			if(nMutation > 0 ){
				for(j=0; j<nMutation; j++){
					pos = rand()%(2*N*nNtrlLoci);
					newPopulation[i].ntrlLoci[pos] =  gsl_rng_uniform_int(r, MAX_NUMBER_OF_INITIAL_NTRL_ALLELES-1) + 1; // ntrl alleles in [1, MAX_NUMBER_OF_INITIAL_NTRL_ALLELES[
				}
			}
			nMutation = 0;
			nMutation = gsl_ran_binomial(r, quantiMutation, 2*N*nQuantiLoci);
			if(nMutation > 0){
				for(j=0; j<nMutation; j++){
					pos = rand()%(2*N*nQuantiLoci);
					currentAllelicEffect = newPopulation[i].quantiLoci[pos];
					minEffect = (1 - RANGE) * currentAllelicEffect;
					if(minEffect<0){
						minEffect = 0.0;
					}
					maxEffect = (1 + RANGE) * currentAllelicEffect;
					if(maxEffect>1.0/2/nQuantiLoci){
					maxEffect = 1.0/2/nQuantiLoci;
					}
					newPopulation[i].quantiLoci[pos] = gsl_ran_flat(r, minEffect, maxEffect);
				}
			}

			pos = 0;
			for(j=0; j<N; j++){	// loop along the individuals to calculate the new femaleAllocation and the number of offsprings
				tmp = 0;
				do{
					newPopulation[i].femaleAllocation[j] += newPopulation[i].quantiLoci[pos];	// sum of allelic effects within individuals
					pos +=1;
					tmp +=1;
				}while(tmp<2*nQuantiLoci);

				newPopulation[i].maleAllocation[j] = 1 - newPopulation[i].femaleAllocation[j];	// set the male allocation as (1 - femaleAllocation)
				newPopulation[i].nOffsprings[j] = floor(fecundity * newPopulation[i].femaleAllocation[j]) + gsl_ran_binomial(r, (fecundity * newPopulation[i].femaleAllocation[j]) - floor(fecundity * newPopulation[i].femaleAllocation[j]), 1);	// set the #of_offsprings produced for each individual

				if(newPopulation[i].sex[j] == 0){ // if heterogametic individual
//				if(newPopulation[i].sexChro[2*j] != newPopulation[i].sexChro[2*j+1]){ // if heterogametic individual
					if(sexualSystem == 1){  // if XY system
						newPopulation[i].femaleAllocation[j] = 0;
						newPopulation[i].maleAllocation[j] = sexAvantage;
						newPopulation[i].nOffsprings[j] = 0;
					}
					if(sexualSystem == 2){  // if ZW system
						newPopulation[i].femaleAllocation[j] = sexAvantage;
						newPopulation[i].maleAllocation[j] = 0;
						double nBabies = fecundity * sexAvantage * 1.0;
						newPopulation[i].nOffsprings[j] = floor(nBabies) + gsl_ran_binomial(r, nBabies - floor(nBabies), 1);
					}
				}
			}
			free(parentalIndexes);
			free(mothers);
			free(fathers);
	}	// end of loop over the nDemes
}

void replacement(gsl_rng* r, Deme* population, Deme* newPopulation, const int nDemes, const int maxIndPerDem, const int nNtrlLoci, const int nQuantiLoci, const int nImmigrants[], int nProducedSeeds[], int extinctionStatus[], const int recolonization, const int generation, const int sexualSystem){
	int i = 0;
	int j = 0;
	int k = 0;
	int compteur = 0;
	int taille = 0;
	int demeTMP = 0;
	int indTMP = 0;
	int posDonneur = 0;	// first position in the emigrant population
	int posReceveur = 0;	// first position in the imigrant population
	int nExtinctedDemes = 0;	// number of extincted demes.
	double* indexOfDemes = NULL;
	double* nProducedSeedsDouble = NULL;

	// count the number of extincted demes
	for(i=0; i<nDemes; i++){
		if(extinctionStatus[i] == 1){
			nExtinctedDemes += 1;
		}
	}

	if(nExtinctedDemes == nDemes){
		printf("Demes are dead at generation: %d\n", generation);
		exit(1);
	}

	indexOfDemes = malloc((nDemes - nExtinctedDemes) * sizeof(double)); // there are (nDemes - nExtinctedDemes) allow to receive or send migrants/recolonizer
	nProducedSeedsDouble = malloc((nDemes - nExtinctedDemes) * sizeof(double));

	if(indexOfDemes == NULL || nProducedSeedsDouble == NULL){
		exit(0);
	}

	for(i=0; i<nDemes; i++){
		if(extinctionStatus[i] != 1){
			indexOfDemes[compteur] = i;
			nProducedSeedsDouble[compteur] = (double)nProducedSeeds[i];
			compteur += 1;
		}
	}

	for(i=0; i<nDemes; i++){
		taille = newPopulation[i].nIndividus + nImmigrants[i];
		if(extinctionStatus[i] == 1){
			taille = recolonization;
		}
		if(taille > maxIndPerDem){
			taille = maxIndPerDem;
		}

		population[i].nIndividus = taille;
		population[i].ntrlLoci = malloc(2 * taille * nNtrlLoci * sizeof(long));
		population[i].quantiLoci = malloc(2 * taille * nQuantiLoci * sizeof(long));
		population[i].sexChro = malloc(2 * taille * sizeof(int));
		population[i].sex = malloc(taille * sizeof(int));
		population[i].femaleAllocation = malloc(taille * sizeof(long));
		population[i].maleAllocation = malloc(taille * sizeof(long));
		population[i].nOffsprings = malloc(taille * sizeof(int));

		if(population[i].ntrlLoci == NULL || population[i].quantiLoci == NULL || population[i].sexChro == NULL || population[i].sex == NULL || population[i].femaleAllocation == NULL || population[i].maleAllocation == NULL || population[i].nOffsprings == NULL){
			exit(0);
		}

		// copy paste newPopulation into population of non extincted demes
		if(extinctionStatus[i] == 0){
			for(j=0; j<(2 * (population[i].nIndividus - nImmigrants[i]) * nNtrlLoci); j++){
				population[i].ntrlLoci[j] = newPopulation[i].ntrlLoci[j];
			}

			for(j=0; j<(2 * (population[i].nIndividus - nImmigrants[i]) * nQuantiLoci); j++){
				population[i].quantiLoci[j] = newPopulation[i].quantiLoci[j];
			}

			for(j=0; j<(2 * (population[i].nIndividus - nImmigrants[i])); j++){
                		population[i].sexChro[j] = newPopulation[i].sexChro[j];
			}

			for(j=0; j<(population[i].nIndividus - nImmigrants[i]); j++){
				population[i].femaleAllocation[j] = newPopulation[i].femaleAllocation[j];
				population[i].maleAllocation[j] = newPopulation[i].maleAllocation[j];
				population[i].nOffsprings[j] = newPopulation[i].nOffsprings[j];
				population[i].sex[j] = newPopulation[i].sex[j];
			}
			// add migrants
			if(nImmigrants[i] > 0){
				double* emigrantDemes = NULL;
				emigrantDemes = malloc(nImmigrants[i] * sizeof(double)); // indexes of demes where immigrants come from
				if(emigrantDemes == NULL){
					exit(0);
				}

				weightedSample(r, indexOfDemes, nProducedSeedsDouble, emigrantDemes, nDemes-nExtinctedDemes, nImmigrants[i]); // return the indexes of nImmigrants[i] emigrant demes from nDemes

				for(j=0; j<nImmigrants[i]; j++){	// loop over the nImmigrants for the deme "i"
					demeTMP = emigrantDemes[j];	// get the emigrant deme
					indTMP = gsl_ran_flat(r, 0, newPopulation[demeTMP].nIndividus);	// randomly choose the emigrant individual from the emigrant deme

					// neutral loci
					posDonneur = 2 * nNtrlLoci * indTMP + 0;
					posReceveur = (2* nNtrlLoci * population[i].nIndividus) - (2 * nNtrlLoci * nImmigrants[i]) + 2 * nNtrlLoci * j;

					for(k=0; k < 2 * nNtrlLoci; k++){	// loop over positions to bring migrants through copy-pasting
						population[i].ntrlLoci[posReceveur + k] = newPopulation[demeTMP].ntrlLoci[posDonneur + k];
					}

					// quantitative loci
					posDonneur = 2 * nQuantiLoci * indTMP + 0;
					posReceveur = (2*nQuantiLoci * population[i].nIndividus) - (2 * nQuantiLoci * nImmigrants[i]) + 2 * nQuantiLoci * j ;

					for(k=0; k < 2 * nQuantiLoci; k++){
						population[i].quantiLoci[posReceveur + k] = newPopulation[demeTMP].quantiLoci[posDonneur + k];
					}

					// sex chromosomes
					posDonneur = 2 * indTMP + 0;
					posReceveur = (2 * population[i].nIndividus) - (2 * nImmigrants[i]) + 2 * j;

					for(k=0; k<2; k++){
						population[i].sexChro[posReceveur + k] = newPopulation[demeTMP].sexChro[posDonneur +k];
					}

					// copy paste femaleAllocation, maleAllocation and nOffsprings
					posDonneur = indTMP;
					posReceveur = population[i].nIndividus - nImmigrants[i] + j;

					population[i].femaleAllocation[posReceveur] = newPopulation[demeTMP].femaleAllocation[posDonneur];
					population[i].maleAllocation[posReceveur] = newPopulation[demeTMP].maleAllocation[posDonneur];
					population[i].nOffsprings[posReceveur] = newPopulation[demeTMP].nOffsprings[posDonneur];

					population[i].sex[posReceveur] = newPopulation[demeTMP].sex[posDonneur];

				}	// end of loop over the n immigrants to the deme "i"

				free(emigrantDemes);
			}	// end of migrant traitement
		}	// end of treatment of non-extincted demes
		if(extinctionStatus[i] == 1){
			double* emigrantDemes = NULL;
			emigrantDemes = malloc(population[i].nIndividus * sizeof(double));
			if(emigrantDemes == NULL){
				exit(0);
			}

			weightedSample(r, indexOfDemes, nProducedSeedsDouble, emigrantDemes, (nDemes - nExtinctedDemes), population[i].nIndividus); // return the indexes of nImmigrants[i] emigrant demes from nDemes. When extinctionStatus==1, all individuals making the deme are imigrants
			for(j=0; j<population[i].nIndividus; j++){	// loop over the individuals to put into extincted demes
				demeTMP = emigrantDemes[j];

				// if we only have hermaphrodites (sexualSystem == 0): all individuals have the same probability of being colonizers 
				if(sexualSystem == 0){
					indTMP = gsl_ran_flat(r, 0, newPopulation[demeTMP].nIndividus);	// randomly choose the colonizer individual from the emigrant deme
				}else{
					// choose the colonizer individual by avoiding unisexuals. Only cosexuals can recolonize 
					double* indexOfIndividuals = NULL; // indexOfIndividuals = [0, 1, 2, ..., N-1] if there are N individuals in the deme i
					double* recolonizer = NULL; // vector of size 1 containing the sampled individual contributing to recolonization
					double* sexDouble = NULL; // convert the vector of (int) sexes [0, 0, 1, 0, 1, ... ] in doubles [0.0, 0.0, 1.0, 1.0, ... ]
					indexOfIndividuals = malloc(newPopulation[demeTMP].nIndividus * sizeof(double)); 
					recolonizer = malloc(1 * sizeof(double));
					sexDouble = malloc(newPopulation[demeTMP].nIndividus * sizeof(double));
					if(indexOfIndividuals == NULL || recolonizer == NULL){
						exit(0);
					}
					for(k=0; k<newPopulation[demeTMP].nIndividus; k++){
						indexOfIndividuals[k] = (double)k;
						sexDouble[k] = (double)newPopulation[demeTMP].sex[k];
					}
				
					weightedSample(r, indexOfIndividuals, sexDouble, recolonizer, newPopulation[demeTMP].nIndividus, 1);
					indTMP = recolonizer[0];
					free(indexOfIndividuals);
					free(recolonizer);
					free(sexDouble);
				}

				// neutral loci
				posDonneur = 2 * nNtrlLoci * indTMP + 0;
				for(k=0; k< 2 * nNtrlLoci; k++){
					population[i].ntrlLoci[2 * nNtrlLoci * j + k] = newPopulation[demeTMP].ntrlLoci[2 * nNtrlLoci * indTMP + k];
				}

				// quantitative loci
				posDonneur = 2 * nQuantiLoci * indTMP + 0;
				for(k=0; k< 2 * nQuantiLoci; k++){
					population[i].quantiLoci[2 * nQuantiLoci * j + k] = newPopulation[demeTMP].quantiLoci[2 * nQuantiLoci * indTMP + k];
				}

		                // sex chromosomes
				posDonneur = 2 * indTMP + 0;
				for(k=0; k< 2; k++){
					population[i].sexChro[2 * j + k] = newPopulation[demeTMP].sexChro[2 * indTMP + k];
				}


				// copy paste femaleAllocation, maleAllocation and nOffsprings
				population[i].femaleAllocation[j] = newPopulation[demeTMP].femaleAllocation[indTMP];
				population[i].maleAllocation[j] = newPopulation[demeTMP].maleAllocation[indTMP];
				population[i].nOffsprings[j] = newPopulation[demeTMP].nOffsprings[indTMP];

				population[i].sex[j] = newPopulation[demeTMP].sex[indTMP];

			}	// end of loop over the individuals to put into extincted demes
			free(emigrantDemes);
		}	// end of treatment of extincted demes
	}	// end of loop over demes

	free(indexOfDemes);
	free(nProducedSeedsDouble);
}

void weightedSample(gsl_rng* r, const double* liste, const double* weights, double* target, const int sizeOfListe, const int nTrials){
	// function that fills the vector 'target' of size 'nTrials' containing the weighted-sampled 'sizeOfListe' elements of the 'liste':
	// weightedSample(gsl_rng* r, {2, 4, 6, 8, 10}, {1.2, 0.6, 0.3, 0.15, 0.05}, target, 5, 20)
	// target = {6, 6, 6, 2, 2, 2, 2, 2, 8, 2, 8, 6, 4, 2, 2, 4, 4, 4, 4, 2}
	// but can also be used for boolean sampling (pile ou face) using:
	// weightedSample(gsl_rng* r, {0, 1}, {1, 1}, target, 2, 1)
	int i = 0;
	unsigned int* n = NULL;
	int* sampledListe = NULL;
	n = malloc(sizeOfListe * sizeof(double));	// will contain the number of succes after K nTrials for each of the sizeOfListe elements of liste
	sampledListe = malloc(nTrials * sizeof(int));	// if n={0, 3, 1, 1}, sampledListe={1, 1, 1, 4, 5}

	gsl_ran_multinomial(r, sizeOfListe, nTrials, weights, n);	// return in 'n' the number of success for the sizeOfListe elements of liste

	int nValues = 0;
	int tmp = 0;
	int tmp2 = 0;
	for(i=0; i<sizeOfListe; i++){ // loop along the list called 'n' resulting from gsl_ran_multinomial
		nValues = n[i];
		if(nValues != 0){
			tmp2 = 0;
			do{
				sampledListe[tmp] = i;
				tmp++;
				tmp2++;
			}while(tmp2 < nValues);
		}
	}

	// shuffle values of the sampledListe
	gsl_permutation* p = gsl_permutation_alloc (nTrials);
	gsl_permutation_init (p);
	gsl_ran_shuffle(r, p -> data, nTrials, sizeof(size_t));

	tmp = 0;
	for(i=0; i<nTrials; i++){
		tmp=gsl_permutation_get(p, i);
		target[i]=liste[sampledListe[tmp]];
	}
	gsl_permutation_free(p);
	free(n);
	free(sampledListe);
}

void libererMemoirePopulation(Deme* population, const int nDemes){
	// free memory taken by the population at the end of each generation.
	int i = 0;
	for(i=0; i<nDemes; i++){
		free(population[i].ntrlLoci);
		free(population[i].quantiLoci);
		free(population[i].sexChro);
		free(population[i].sex);
		free(population[i].femaleAllocation);
		free(population[i].maleAllocation);
		free(population[i].nOffsprings);
	}
}

void writeNindividuals(const Deme* population, const int nDemes, const double extinction, const double migration, const int seed){
	int i = 0;

	char nameOfFileNInd[100];
	char nameOfFileFemAlloc[100];
	sprintf(nameOfFileNInd, "indOverTime_e%f_i%f_seed%d.txt", extinction, migration, seed);
	sprintf(nameOfFileFemAlloc, "femAllocOverTime_e%f_i%f_seed%d.txt", extinction, migration, seed);

	FILE* fichierNInd = NULL;
	FILE* fichierFemAlloc = NULL;
	fichierNInd =  fopen(nameOfFileNInd, "a");
	fichierFemAlloc =  fopen(nameOfFileFemAlloc, "a");
	if(fichierNInd != NULL && fichierFemAlloc != NULL){
		for(i=0; i<nDemes; i++){
			fprintf(fichierNInd, "%d ", population[i].nIndividus);
			fprintf(fichierFemAlloc, "%f ", gsl_stats_mean(population[i].femaleAllocation, 1, population[i].nIndividus));
		}
		fprintf(fichierNInd, "\n");
		fclose(fichierNInd);
		fprintf(fichierFemAlloc, "\n");
		fclose(fichierFemAlloc);
	}
}

void afficherPopulation(Deme* population, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const int sexualSystem){
	// called to print some informations about population in a debug mode
	int i = 0;
	int j = 0;
	int k = 0;
	int sexGenotype = 0;
	char sex = 'H';
	for(i=0; i<nDemes; i++){
		for(j=0; j<population[i].nIndividus; j++){
			if(population[i].sexChro[2*j] == population[i].sexChro[2*j + 1]){
				sexGenotype = 0;	// homogametic
			}
			if(population[i].sexChro[2*j] != population[i].sexChro[2*j + 1]){
				sexGenotype = 1;	// heterogametic
			}
			if(sexualSystem == 0){
				sex = 'H';
			}
			if(sexualSystem == 1){
				if(sexGenotype == 0){
					sex = 'H';
				}
				if(sexGenotype == 1){
					sex = 'M';
				}
			}
			if(sexualSystem == 2){
				if(sexGenotype == 0){
					sex = 'H';
				}
				if(sexGenotype == 1){
					sex = 'F';
				}
			}
			
			printf("Deme: %d Ind: %d Ntrl: ", i, j);
			for(k=0; k<(2*nNtrlLoci); k++){
				printf("%ld ", population[i].ntrlLoci[2*j*nNtrlLoci+k]);
			}
			printf(" Quanti: ");
			for(k=0; k<(2*nQuantiLoci); k++){
				printf("%.4lf ", population[i].quantiLoci[2*j*nQuantiLoci+k]);
			}
		printf("femAlloc: %.4lf nOffs: %d sex: %d sex2: %c\n", population[i].femaleAllocation[j], population[i].nOffsprings[j], population[i].sex[j], sex);
		}
	}
}

void genePop(Deme* population, const int nDemes, const int nNtrlLoci, const int seed, int time){
	//	generates an input for genepop (Rousset), launch it and clean the output to save space on the cluster
	int i = 0;
	int j = 0;
	int k = 0;
	int allele = 0;
	int cntLoci = 0;
	char nameOfGenePopFile[100];
	char nameOfSettingFile[100];
	char nameOfROutputFile[100];
//	char commandLineOne[100]; // launch genepop
//	char commandLineTwo[200]; // treat the genepop's output
	char commandLineThree[200]; // clean the tmp files
	char commandLineDiveRsity[200]; // R command calling diveRsity 
	sprintf(nameOfGenePopFile, "genepop_%d_%d.txt", time, seed);
	sprintf(nameOfSettingFile, "setting_%d_%d.txt", time, seed);
	sprintf(nameOfROutputFile, "output_diveRsity_%d_%d", time, seed);

//	sprintf(commandLineOne, "Genepop settingsFile=%s Mode=Batch >/dev/null", nameOfSettingFile); // external call of genepop (Rousset)
//	sprintf(commandLineTwo, "tail -n%d %s.FST | grep 'Locus'>tmp_%d.txt; tail -n%d %s.FST | grep 'All' >>tmp_%d.txt; mv tmp_%d.txt %s.FST", nNtrlLoci + 5, nameOfGenePopFile, seed, nNtrlLoci + 5, nameOfGenePopFile, seed, seed, nameOfGenePopFile); // external formating of genepop's output.
	sprintf(commandLineThree, "rm -rf cmdline.txt fichier.in %s %s", nameOfGenePopFile, nameOfSettingFile);

	sprintf(commandLineDiveRsity, "diveRsity.R input=%s output=%s", nameOfGenePopFile, nameOfROutputFile);

	FILE* fichierGenePop = NULL;
	FILE* settingFile = NULL;

	fichierGenePop = fopen(nameOfGenePopFile, "a");	// creating genepop file with data
	if(fichierGenePop != NULL){	// writing in file

		fprintf(fichierGenePop, "Simulated data\n");
	
		for(i=0; i<nNtrlLoci; i++){
			fprintf(fichierGenePop, "Locus%d\n", i);
		}

		for(i=0; i<nDemes; i++){	// loop over demes: start
			fprintf(fichierGenePop, "Pop\n");
			for(j=0; j<population[i].nIndividus; j++){	// loop over individuals of deme 'i': start
				fprintf(fichierGenePop, "Ind%d, ", j);
				for(k=0; k<(2*nNtrlLoci); k++){
					cntLoci += 1;
					allele = population[i].ntrlLoci[2*j*nNtrlLoci+k];
					if(allele < 10){
						fprintf(fichierGenePop, "00%d", allele);
					}
					if(allele >= 10 && allele < 100){
						fprintf(fichierGenePop, "0%d", allele);
					}
					if(allele >= 100){
						fprintf(fichierGenePop, "%d", allele);
					}
					if(cntLoci == 2){
						cntLoci = 0;
						fprintf(fichierGenePop, " ");
					}
				}
				fprintf(fichierGenePop, "\n");
			}	// loop over individuals of deme 'i': end
		}	// loop over demes: end
	}	// end of writing in file
	fclose(fichierGenePop);	// input file for genepop is generated

	settingFile = fopen(nameOfSettingFile, "a");	// creating setting file for genepop
	if(settingFile != NULL){
		fprintf(settingFile, "GenepopInputFile=%s\nMenuOptions=6.1\n", nameOfGenePopFile);
	}
	fclose(settingFile);

//	int testCommandLineOne = system(commandLineOne); // call genepop
//	if(testCommandLineOne == -1){
//		exit(1); // check the returned value
//	}

	int testNameOfROutputFile = system(commandLineDiveRsity); // call diveRsity R
	if(testNameOfROutputFile == -1){
		exit(1); // check the returned value
	}

//	int testCommandLineTwo = system(commandLineTwo); // reformat the genepop's output
//	if(testCommandLineTwo == -1){
//		exit(1); // check the returned value
//	}
	
	int testCommandLineThree = system(commandLineThree); // remove genePop formated file and setting file
	if(testCommandLineThree == -1){
		exit(1);
	}

}

void statisticsPopulations(Deme* population, const int nDemes, const int maxIndPerDem, const int nQuantiLoci, const int fecundity, const double migration, const double extinction, const int recolonization, const int sexualSystem, const double sexAvantage, const int seed, int time){
	// function that calculates the mean female allocation, its standard deviation and the percentage of cosexuals in the metapopulation
	int i = 0;
	int j = 0;
	int cnt = 0;
	int cnt2 = 0;
	double fstValue = 0.0;
	double fstValueDensity = 0.0;
	double meanAllocFemale = 0.0;
	double sdAllocFemale = 0.0;
	double meanAllocFemaleCosexual = 0.0;
	double sdAllocFemaleCosexual = 0.0;
	double cosexualProportion = 0.0;
	int nIndividusTotal = 0;

	char nomFichierSortie[200];
	sprintf(nomFichierSortie, "output_%d.txt", seed);
	FILE* fichierSortie = NULL;

	fichierSortie = fopen(nomFichierSortie, "r");
	if(fichierSortie == NULL){
		fichierSortie = fopen(nomFichierSortie, "a");
		fprintf(fichierSortie, "nDemes\tnIndMaxPerDeme\tNtot\tnQuantiLoci\tfecundity\tmigRate\textRate\trecolonization\tatGeneration\tsexSystem\tsexAvantage\tseed\tmeanFemAlloc\tsdFemAlloc\tmeanFemAllocCosexual\tsdFemAllocCosexual\tcosexualProportion\texpFST_Nmax\texpFST_Nobs\n");
		fclose(fichierSortie);
	}else{
		fclose(fichierSortie);
	}
	
	fichierSortie = fopen(nomFichierSortie, "a");

	if(fichierSortie != NULL){	
	
		for(i=0; i<nDemes; i++){
			nIndividusTotal += population[i].nIndividus;
		}	
	
		double* allocFemale = NULL; // female allocation in the whole metapopulation
		double* allocFemaleCosexual = NULL; // female allocation of cosexuals only. allocFemale = allocFemaleCosexual if sexualSystem == 0

		for(i=0; i<nDemes; i++){
			for(j=0; j<population[i].nIndividus; j++){
				cosexualProportion += population[i].sex[j]; // sex[j] = 0 if unisexual; sex[j] + 1 if cosexual 
			}
		}

		allocFemale = malloc(nIndividusTotal * sizeof(double));
		allocFemaleCosexual = malloc(cosexualProportion * sizeof(double));

		for(i=0; i<nDemes; i++){
			for(j=0; j<population[i].nIndividus; j++){
				allocFemale[cnt] = population[i].femaleAllocation[j];
				if(population[i].sex[j] == 1){
					allocFemaleCosexual[cnt2] = population[i].femaleAllocation[j];
					cnt2 += 1;
				}
				cnt += 1;
			}
		}

		meanAllocFemale = gsl_stats_mean(allocFemale, 1, nIndividusTotal);
		sdAllocFemale = gsl_stats_sd(allocFemale, 1, nIndividusTotal);
		meanAllocFemaleCosexual = gsl_stats_mean(allocFemaleCosexual, 1, (int) cosexualProportion);
		sdAllocFemaleCosexual = gsl_stats_sd(allocFemaleCosexual, 1, (int) cosexualProportion);
		
		cosexualProportion = cosexualProportion / nIndividusTotal;
	
		fstValue = fst(maxIndPerDem, extinction, recolonization, migration);
		fstValueDensity = fst((int) nIndividusTotal/(1.0*nDemes), extinction, recolonization, migration);
	
		free(allocFemale);
		free(allocFemaleCosexual);
	
		fprintf(fichierSortie, "%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%d\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", nDemes, maxIndPerDem, nIndividusTotal, nQuantiLoci, fecundity, migration, extinction, recolonization, time, sexualSystem, sexAvantage, seed, meanAllocFemale, sdAllocFemale, meanAllocFemaleCosexual, sdAllocFemaleCosexual, cosexualProportion, fstValue, fstValueDensity);
		fclose(fichierSortie);
	}
}

void checkCommandLine(int argc){
	if(argc != 15){
		printf("\n%sThe number of provided arguments is not correct.%s\nYou provided %s%d%s argument while %s14%s are expected:\n\t\
		%s1.%s  Number of demes (>0)\n\t\
		%s2.%s  Max number of individuals per deme (>0)\n\t\
		%s3.%s  Number of generations (>0)\n\n\t\
		%s4.%s  Number of neutral loci (>=0)\n\t\
		%s5.%s  Neutral mutation rate (in [0-1]))\n\t\
		%s6.%s  Number of quantitative loci (>0)\n\t\
		%s7.%s  Quantitative mutation rate (in [0-1])\n\n\t\
		%s8.%s  Max number of offsprings per hermaphrodite (>0)\n\n\t\
		%s9.%s  Immigration rate (Poisson distributed; >=0)\n\t\
		%s10.%s Extinction rate (Binomialy distributed; in [0-1])\n\t\
		%s11.%s Number of individuals recolonizing an extincted deme (>0)\n\n\t\
		%s12.%s sexualSystem is equal to 0 if autosomal, equal to 1 if XY and equal to 2 if ZW\n\t\
		%s13.%s Sexual effects of heterogametic sex (if equal to 1.5 in XY system, males have a 50 percent advantage to sire available ovules). Required but neglected if sexualSystem == 0\n\t\
		%s14.%s Seed for the random generator (>0)\n\n\
		%s\tExample:%s ./quantiSex 100 100 100 10 0.0001 1 0.00001 100 1 0.1 1 0 2 123\n\n\
		version: %s\n\n\t\tdependencies: \t%s\n\n", KRED, STOP, KRED, argc-1, STOP, KRED, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KRED, STOP, VERSION, DEPENDENCY);

		exit(0);
	}
}

double fst(const int maxIndPerDem, const double extinction, const int recolonization, const double migration){
	double res = 0.0;
	double numQ = 0.0;
	double denomQ = 0.0;
//	double phi = 0.0;
	double qr = 0.0;
	double migrationProportion = migration / maxIndPerDem;
	
	numQ = 1/(2.0 * maxIndPerDem) + extinction/(2.0 * recolonization) - extinction/(2.0 * recolonization * 2.0 * maxIndPerDem);
	denomQ = 1 - (1 - 1/(2.0 * maxIndPerDem)) * ((1 - migrationProportion)*(1 - migrationProportion) * (1 - extinction) + extinction * (1 - 1/(2.0 * recolonization)) * 1/(2.0 * recolonization -1));
//	phi = 1/(2.0 * recolonization -1);
	qr = numQ/denomQ;
	res = (qr - 1/(2.0 * maxIndPerDem)) * (2 * maxIndPerDem)/(2.0 * maxIndPerDem -1);

	return(res);
}

/*# diveRsityR
#!/usr/bin/env Rscript
./diveRsity.R input=nameOfGenePopFile output=nameOfROutputFile
library(diveRsity)
options(warn=-1)
for(i in commandArgs()){
	tmp = strsplit(i, "=")
	if(tmp[[1]][1] == "input"){input = tmp[[1]][2]}
	if(tmp[[1]][1] == "output"){output = tmp[[1]][2]}
}

a=diffCalc(input, fst=T, pairwise=F, outfile=output)

*/

