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
//	int *tableau;
	int *Random;
	int *Db;
	int *Dmax;
	int* RDRk;
	int tailleTableau, Wn, tailleRDR;
	int pow, borne, seed;
	int nbreTest, i;
	mpz_t k, powWn;
	long int t;

    if ( argc != 5) 
     { 
       printf( "Syntax: RDR <nombre de Tests> <taille de D> <taille de k> <borne de D> \a\n " ); 
       exit(1); 
     } 


	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, time(NULL));

					//////////////////
					/* Méthode RDR */
					//////////////////
	t=cputime();
	tailleTableau = strtol(argv[2],NULL,10); // nombre de coeff dans D
	seed = 5; // graine du générateur pseudo-aléatoire
	nbreTest = strtol(argv[1],NULL,10); // nombre de test que l'on effectue 
	int n = strtol(argv[3],NULL,10); // taille de k
	borne = strtol(argv[4],NULL,10); // valeur maximale du plus grand coefficient de D


/////////////////////////////////////////////////////////////////////////////////////////
/* Création du tableau D avec nos coefficients (attention 1 doit être dans le tableau) */
/////////////////////////////////////////////////////////////////////////////////////////
/*
	tableau = malloc(tailleTableau * sizeof(int));
	initialiserTableau(tableau,tailleTableau);

	tableau[0] = 1;
	tableau[1] = 3;
	tableau[2] = 23;
	tableau[3] = 27;
*/
////////////////////////////////////////////////////////
/* Répétition nbreTest fois du processus pour le test */
////////////////////////////////////////////////////////

	for (i=0;i<nbreTest;i++){

//	printf("\n--------------------------\nTest n°%d :\n---------------------------\n",i); 

///////////////////////////////////////////////////////
/* Définition de k et de powWn (sous forme de mpz_t) */
///////////////////////////////////////////////////////

	mpz_init(k);
	mpz_urandomb(k,state,n); 
//	mpz_init_set_ui(k,706440);

//	gmp_printf("\nk is : %Zd\n",k);	

	mpz_init(powWn);

///////////////////////////////////////////////////////////
/* Création d'un tableau random qui fera notre tableau D */
///////////////////////////////////////////////////////////

	Random = malloc(tailleTableau * sizeof(int));
	generateRandomD(Random,tailleTableau,borne,seed+i);
	afficherTableau(Random,tailleTableau);
//	Random[0] = 1;
//	Random[1] = 3;
//	Random[2] = 23;
//	Random[3] = 27;

//////////////////////////////////////////////////////////////////////////
/* Calcule de Wn puis création du fameux tableau Dmax de taille 2^(Wn+2)*/
//////////////////////////////////////////////////////////////////////////

	Wn = calculWn(Random,tailleTableau);

	pow = pow2i(Wn+2);

	mpz_set_ui(powWn,pow); 

	Dmax = malloc(pow * sizeof(int));
	initialiserTableau(Dmax,pow);

////////////////////////////////
/* Création du tableau Dbarre */
////////////////////////////////

	Db = malloc(2*tailleTableau * sizeof(int));
	initialiserTableau(Db,tailleTableau*2);

	remplirDb(Db,Random,tailleTableau,Wn);

/////////////////////////
/* Remplissage de Dmax */
/////////////////////////

	remplirDmaxAvecDb(Dmax,Db,Random,tailleTableau);
	remplirDmax(Dmax,Db,tailleTableau,Wn,pow);

////////////////////////
/* Calcul de digit(k) */
////////////////////////

//	digitDk = digitD(Dmax,k,powWn);

///////////////////////////////////
/* Calcul de la RDR et affichage */
///////////////////////////////////

	tailleRDR = mpz_sizeinbase(k,2)*2;
	RDRk = malloc(tailleRDR * sizeof(int));
	RDRk = RDR(Dmax,k,Wn);

/////////////////////////
/* Affichage et autres */
/////////////////////////

	printf("Wn is %d\n",Wn);

//	afficherTableau(Db,2*tailleTableau);
//	afficherTableau(Dmax,pow);
	affichageRDR(RDRk,k);

	free(Random);
	free(Dmax);
	free(Db);
	free(RDRk);

	}

	gmp_randclear(state);

	printf("\nNombre de tests : %d\n",nbreTest);
	printf("\nTaille de D : %d\n",tailleTableau);
	printf("\nTaille de k : %d\n",n);
	printf("\n\nRDR time: %ld\n",cputime()-t);
	printf("\n-------------------------------\n");

	return 0;
}
