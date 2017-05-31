#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <inttypes.h>
#include <time.h>
#include <gmp.h>
#include <unistd.h>
#include <sys/time.h>
#include "RDRgmp.h"

long int cputime (){
  struct timeval time;
  struct timezone z;
  
  gettimeofday (&time, &z);
  return (long int) time.tv_sec * 1000 + time.tv_usec / 1000;
}

int main(int argc, char * argv[]){
//	int* RDwNAF;
	int w, tailleRDwNAF, tailleDmax;
	int nbreTest, i;
	mpz_t k;
	long int t;

    if ( argc != 4) 
     { 
       printf( "Syntax: wNAF <nombre de Tests> <taille de D> <taille de k> \a\n " ); 
       exit(1); 
     } 

					//////////////////
					/* Méthode wNAF */
					//////////////////

	nbreTest = strtol(argv[1],NULL,10);
	w = strtol(argv[2],NULL,10);
//	w = 64; // taille de la fenetre w
	tailleDmax = 4*w;
//	nbreTest = 10000;
	t=cputime();


	int n = strtol(argv[3],NULL,10);

	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, time(NULL));
//	mp_bitcnt_t n = 256;

////////////////////////////////////////////////////////
/* Répétition nbreTest fois du processus pour le test */
////////////////////////////////////////////////////////

	for (i=0;i<nbreTest;i++){

//	printf("\n--------------------------\nTest n°%d :\n---------------------------\n",i); 

///////////////////////////////////////////////////////
/* Définition de k (sous forme de mpz_t) */
///////////////////////////////////////////////////////

	mpz_init(k);
	mpz_urandomb(k,state,n); 

//	gmp_printf("\nk is : %Zd\n",k);

////////////////////////////
/* Calcul du Dmax de wNAF */
////////////////////////////

	int* DmaxwNAF=calloc(tailleDmax,sizeof(int));

	wNAFDmax(DmaxwNAF,tailleDmax);

/////////////////////////////////////////////////////////
/* Calcul de la représentation wNAF de k sur fenetre w */
/////////////////////////////////////////////////////////

	tailleRDwNAF = mpz_sizeinbase(k,2)*2;
	int* RDwNAF=calloc(tailleRDwNAF,sizeof(int));

	RDwNAF = wNAF(DmaxwNAF, tailleDmax,k);

////////////////////////////////////
/* Affichage de la représentation */
////////////////////////////////////

//	afficherTableau(DmaxwNAF,tailleDmax);
//	affichageRDR(RDwNAF,k);

	free(DmaxwNAF);
	free(RDwNAF);
	}

	gmp_randclear(state);

	printf("\nNombre de tests : %d\n",nbreTest);
	printf("\nTaille de D : %d\n",w);
	printf("\nTaille de k : %d\n",n);
	printf("\n\nwNAF time: %ld\n",cputime()-t);
	printf("\n-------------------------------\n");

	return 0;
}

