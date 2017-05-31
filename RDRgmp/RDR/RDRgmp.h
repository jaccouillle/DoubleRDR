int pow2i(int i); // on calcule 2^i
int tailleInt(int a); // donner la taille de a en binaire
int maxTableau(int* tableau, int taille); // calcule le maximum d'un tableau
int calculWn(int* D, int tailleD); // on calcule le paramètre Wn de la fonction digit
void initialiserTableau(int *tableau, int tailleTableau); // initialise une tableau de mpz_t
void remplirDb(int* Db, int* D, int tailleD, int Wn); // on calcule Dbarre = D U {2^(Wn+2)-d, d in D}
void afficherTableau(int* tableau, int tailleTableau); // on affiche un tableau de mpz_t
void remplirDmaxAvecDb(int* Dmax, int* Db, int* D, int tailleD); // on remplit Dmax avec Db (1ère étape)
int calculd(int Dbj, int powi, int k, int powWn); // on calcule les d de l'algo superDmax
void remplirDmax(int* Dmax, int* Db, int tailleD, int Wn, int tailleDmax); // on remplit Dmax
int appartientTableau(int element, int* tableau, int taille); // renvoie 1 si element est dans le tableau
void generateRandomD(int* Random, int taille, int borne, int seed); // on génère un ensemble D aléatoirement (avec 1 dedans)
int digitD(int* Dmax, mpz_t k, mpz_t powWn);
int digitDparticular(int* Dmax, mpz_t k, int Wn);
int* RDR(int* Dmax, mpz_t k, int Wn);
void affichageRDR(int* RDR, mpz_t k);
void wNAFDmax(int* Dmax, int taille);
int* wNAF(int* Dmax, int tailleDmax, mpz_t k);
