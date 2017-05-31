#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <inttypes.h>
#include <time.h>
#include <gmp.h>
#include "RDRgmp.h"

////////////////////////
/* fonction puissance */
////////////////////////

int pow2i(int i){ // on calcule 2^i
	int p = 1, j;

	for (j=0;j<i;j++)
		p=p*2;

	return p;
}	

/////////////////////////////////////////////////////
/* fonction qui calcule la taille binaire d'un int */
/////////////////////////////////////////////////////

int tailleInt(int a){ // donner la taille de a en binaire
	int i = 1, b = a;
	while (b > 1 && i < 20000){
		b = b/2;
		i++;
	}
	return i;
}

/////////////////////////////////////////////////////////
/* fonction qui calcule le maximum d'un tableau de int */
/////////////////////////////////////////////////////////

int maxTableau(int* tableau, int taille){ // calcule le maximum d'un tableau
	int i, c = 0;

	for (i=0;i<taille;i++){
		if (tableau[i] > c)
			c = tableau[i];
	}

	return c;
}

//////////////////////////////////////////////////
/* fonction qui calcule Wn de la fonction digit */
//////////////////////////////////////////////////

int calculWn(int* D, int taille){ // on calcule le paramètre Wn de la fonction digit
	int max, Wn;

	max = maxTableau(D,taille);
	Wn = tailleInt(max)-1;

	return Wn;
}

///////////////////////////////////////////////
/* fonction qui initialise un tableau de int */
///////////////////////////////////////////////

void initialiserTableau(int* tableau, int tailleTableau){ // initialise une tableau de int
	int i;

	for (i=0;i<tailleTableau;i++)
		tableau[i] = 0;
}

///////////////////////////////////////////////////////////////////////////////
/* fonction qui créer le tableau Dbarre = D U {2^(Wn+2)-d,d appartenant à D} */
///////////////////////////////////////////////////////////////////////////////

void remplirDb(int* Db, int* D, int tailleD, int Wn){ // on calcule Dbarre = D U {2^(Wn+2)-d, d in D}
	int k;
	int pow = 1;

	for (k=0;k<Wn+2;k++)
		pow = pow*2;

	for (k=0;k<tailleD;k++){
		Db[k] = D[k];
		Db[tailleD+k] = pow - D[k];
	}
}

////////////////////////////////////////////
/* fonction qui affiche un tableau de int */
////////////////////////////////////////////

void afficherTableau(int* tableau, int tailleTableau){ // on affiche un tableau de int
	int i;

	printf("\n Les valeurs du tableau sont : \n");
	for (i=0;i<tailleTableau;i++)
		printf("%d,",tableau[i]);
	printf("\n");

}

////////////////////////////////////////////////////////////////////////
/* fonction qui complète Dmax avec D et -D pour les indices de Dbarre */
////////////////////////////////////////////////////////////////////////

void remplirDmaxAvecDb(int* Dmax, int* Db, int* D, int tailleD){ // on remplit Dmax avec Db (1ère étape)
	int i;

	for (i=0;i<tailleD;i++)
		Dmax[Db[i]] = D[i];

	for (i=0;i<tailleD;i++)
		Dmax[Db[tailleD+i]] = -D[i];
}

//////////////////////////////////////////////////////////////////////////
/* fonction qui calcule chaque indice que l'on rajoute dans l'algo Dmax */
//////////////////////////////////////////////////////////////////////////

int calculd(int Dbj, int powi, int k, int powWn){ // on calcule d = Db[j] + 2^i(2k+1) Mod 2^(Wn+2)
	int d = 0, mult, mod;

	mult = powi*(2*k+1);

	mod = mult + Dbj;

	d = mod%powWn;

	return d;
}

/////////////////////////////////////////
/* fonction qui calcule le fameux Dmax */
/////////////////////////////////////////

void remplirDmax(int* Dmax, int* Db, int tailleD, int Wn, int tailleDmax){ // on remplit Dmax
	int i = Wn+1,j,k,borneK, compteurTableau = 0, Dbj, powWn, powi;
	int d;

	powWn = pow2i(Wn+2); // 2^(Wn+2) utile dans tout l'algo pour l'indiçage
	powi = powWn/2; // 2^(Wn+1) utile pour le début i=Wn+1 (d = Db[j] + 2^(Wn+1) pour tout j)
	borneK = 1; // utile pour le début i=Wn+1

	while (i>2 && compteurTableau<tailleDmax){
		j=2*tailleD-1;

		while (j>-1){
			Dbj = Db[j]; // Db[j] reste identique quand j reste identique
			for (k=0;k<borneK;k++){
				d=calculd(Dbj,powi,k,powWn);

				if (Dmax[d]==0){
					compteurTableau=compteurTableau+2;
					if (j<tailleD)
						Dmax[d] = Db[j];
					else
						Dmax[d] = -Db[j-tailleD];
				}
			}
			j--;
		}
		borneK = borneK*2; // pour dire que l'on fera 2 fois plus d'addition qu'a l'étape précédente
		powi = powi/2; // pour dire que l'on ajoutera des nombres 2 fois plus petit qu'à l'étape précédente
		i--;
	}

	if (compteurTableau<tailleDmax){
		for (i=0;i<tailleDmax/4;i++){
			if (Dmax[4*i+1] == 0)
				Dmax[4*i+1] = 1;
			if (Dmax[4*i+3] == 0)
				Dmax[4*i+3] = -1;
		}
	}
}

/////////////////////////////////////////////////////////////////
/* fonction qui renvoie 1 si element est dans tableau, 0 sinon */
/////////////////////////////////////////////////////////////////

int appartientTableau(int element, int* tableau, int taille){ // renvoie 1 si element est dans le tableau
	int j, r = 0;

	for (j=0;j<taille;j++){
		if (tableau[j] == element)
			r = 1;
	}

	return r;
}

////////////////////////////////////////////////////////////////////////////////////////
/* fonction qui génére aléatoirement un ensemble D contenant 1 et d'éléments <= borne */
////////////////////////////////////////////////////////////////////////////////////////

void generateRandomD(int* Random, int taille, int borne, int seed){ // on génère un ensemble D aléatoirement (avec 1 dedans)
	int i,r;

	Random[0] = 1;

	srand(time(0)+seed);
	for (i=1;i<taille;i++){
		r = rand();
		r = (2*(r%borne))%borne + 3;
		if (appartientTableau(r,Random,taille) == 0 && r<borne)
			Random[i] = r;
		else i--;
	}
}	

////////////////////
/* fonction digit */
////////////////////

int digitD(int* Dmax, mpz_t k, mpz_t powWn){
	mpz_t kmod;
	int kmodint, digit = 0, parite;

	mpz_init(kmod);

	parite = mpz_tstbit(k,0);

	if (parite == 1){
		mpz_mod(kmod,k,powWn);
		kmodint = mpz_get_ui(kmod);

		digit = Dmax[kmodint];
	}

	return digit;
}

/////////////////////////////////////////////////
/* Calcul de digit au cas où di = digitD(k) > k*/
/////////////////////////////////////////////////

int digitDparticular(int* Dmax, mpz_t k, int Wn){
	int kint = mpz_get_ui(k);
	int w = 1, i, digit = kint+1;

	for (i=0;i<Wn;i++)
		w = w*2;

	while (digit > kint || -digit > kint){
		kint = kint%w;
		digit = Dmax[kint];
		w = w/2;
	}

	return digit;
}

///////////////////////////
/* Calcul de la RDR de k */
///////////////////////////

int* RDR(int* Dmax, mpz_t k, int Wn){
	mpz_t interk, zero, deux, powWn, subk, divk;
	int cmp = 1, i = 0, digit;
	int j, tailleRDR;
	int taille = mpz_sizeinbase(k,2);
	int* representationInv=calloc(taille*2,sizeof(int));

	mpz_init(interk);
	mpz_set(interk,k);

	mpz_init(subk);
	mpz_init(divk);
	mpz_init(zero);
	mpz_init_set_ui(deux,2);
	mpz_init(powWn);

	mpz_pow_ui(powWn,deux,Wn+2);

	while (cmp != 0 && i<taille*2){
		digit = digitD(Dmax,interk,powWn);

		if (mpz_cmpabs_d(interk,digit) < 0) // au cas ou k < digit
			digit = digitDparticular(Dmax,interk,Wn);

		representationInv[i] = digit;

		if (digit > -1)
			mpz_sub_ui(subk,interk,digit); // calcule subk = (ki-digit)
		else mpz_add_ui(subk,interk,-digit); // calcule subk = (ki+(-digit)) avec -digit > 0

		mpz_fdiv_q(divk,subk,deux); // calcule divk = subk/2 = (ki-digit)/2
		mpz_set(interk,divk);

		cmp = mpz_cmp(interk,zero); // on check si ki différent de 0

		i++;
	}

	tailleRDR = i;

	int* RDRk=(int*)malloc(sizeof(int)*(i+1));

	for (j=0;j<tailleRDR;j++)
		RDRk[j] = representationInv[tailleRDR-1-j];

	RDRk[i] = 666;

	free(representationInv);

	return RDRk;
}

//////////////////////////////
/* Affichage de la RDR de k */
//////////////////////////////

void affichageRDR(int* RDR, mpz_t k){
	int i = 0, j;
	int taille = mpz_sizeinbase(k,2);

	gmp_printf("\nLa RDR de %Zd de taille %d est :\n",k, taille);
	while (i<taille*2 && RDR[i] != 666){
		printf("%d,",RDR[i]);
		i++;
	}
	gmp_printf("\n%Zd = ",k);
	for (j=0;j<i;j++){
		if (RDR[j] != 0)
			printf("%d*2^(%d)+",RDR[j],i-1-j);
	}
}

////////////////////////////////////////////
/* Fonction qui calcule le Dmax pour wNAF */
////////////////////////////////////////////

void wNAFDmax(int* Dmax, int taille){
	int i;

	for (i=0;i<taille;i++){
		if (i%2 != 0)
			Dmax[i] = i;
		else Dmax[i] = 0;
	}
}

///////////////////////////////////////////////////////////
/* Fonction qui effectue une représentation grâce à wNAF */
///////////////////////////////////////////////////////////

int* wNAF(int* Dmax, int tailleDmax, mpz_t k){
	mpz_t zero, deux, interk, subk, divk;
	int cmp = 1, i = 0, ki, taille, pow, digit, tailleRDwNAF, j;

/* On crée le tableau de la représentation wNAF */

	taille = mpz_sizeinbase(k,2);

	int* representationInv=calloc(taille*2,sizeof(int));

	mpz_init(zero);
	mpz_init_set_ui(deux,2);
	mpz_init(interk);
	mpz_set(interk,k);
	mpz_init(subk);
	mpz_init(divk);
	pow = tailleDmax;
	while (cmp != 0 && i<taille*2){
		ki = mpz_fdiv_ui(interk,pow);
		digit = Dmax[ki];

		representationInv[i] = digit;

		mpz_sub_ui(subk,interk,digit);
		mpz_fdiv_q(divk,subk,deux);
		mpz_set(interk,divk);

		cmp = mpz_cmp(interk,zero);
		i++;
	}

	tailleRDwNAF = i;

	int* RDwNAFk=(int*)malloc(sizeof(int)*(tailleRDwNAF+1));

	for (j=0;j<tailleRDwNAF;j++)
		RDwNAFk[j] = representationInv[tailleRDwNAF-1-j];
	RDwNAFk[tailleRDwNAF] = 666;

	free(representationInv);

	return RDwNAFk;
}


