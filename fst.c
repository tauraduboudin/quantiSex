#include <stdio.h>
#include <stdlib.h>

// estimate the expected Fst value in a metapoluation following:
// Effective size in simple metapopulation models ; F Rousset ; Heredity ; 2003
double fst(const int nDemes, const int maxIndPerDem, const double extinction, const int recolonization, const double migration);

int main(int argc, char *argv[]){
	const int nDemes = atoi(argv[1]);
	const int maxIndPerDem = atoi(argv[2]);
	const double extinction = atof(argv[3]);
	const int recolonization = atoi(argv[4]);
	const double migration = atof(argv[5]);

	const double fstValue = fst(nDemes, maxIndPerDem, extinction, recolonization, migration);

	printf("fst = %lf\n", fstValue);
	return(0);
}

double fst(const int nDemes, const int maxIndPerDem, const double extinction, const int recolonization, const double migration){
	double res = 0.0;
	double numQ = 0.0;
	double denomQ = 0.0;
	double phi = 0.0;
	double qr = 0.0;
	double migrationProportion = migration / maxIndPerDem;
	
	numQ = 1/(2.0 * maxIndPerDem) + extinction/(2.0 * recolonization) - extinction/(2.0 * recolonization * 2.0 * maxIndPerDem);
	denomQ = 1 - (1 - 1/(2.0 * maxIndPerDem)) * ((1 - migrationProportion)*(1 - migrationProportion) * (1 - extinction) + extinction * (1 - 1/(2.0 * recolonization)) * 1/(2.0 * recolonization -1));
	phi = 1/(2.0 * recolonization -1);
	qr = numQ/denomQ;
	res = (qr - 1/(2.0 * maxIndPerDem)) * (2 * maxIndPerDem)/(2.0 * maxIndPerDem -1);

	return(res);
}

