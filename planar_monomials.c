/* Compile with: g++ -O2 -pthread -march=native -funroll-loops planar_monomials.c -o planar_monomials -lntl -lgmp -lm
Author: C. Beierle, P. Felke -- Nov 2021

This program checks the planarity of monomials over a finite field
*/

#include <time.h>
#include <math.h>
#include <algorithm>
#include <set>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>

//if this is set, the code uses the classification of planar monomials over the finite fields F_{p^i}, i=5,...,8
#define USE_DIVISORS_5_TO_8

using namespace std;
using namespace NTL;

unsigned long N_SAMPLES;

unsigned long P;
unsigned long N;
unsigned long Q;
unsigned long P2;
unsigned long P3;
unsigned long P4;

unsigned long P5;
unsigned long P6;
unsigned long P7;
unsigned long P8;

int runtime;

// returns 1 if a is the smallest value in the set {p^i*a mod p^n-1, i=0,...,n-1}
int is_smallest_in_equivalence_class(unsigned long a) {
  unsigned long b;
  for (int i=1; i<N; i++) {
    b = (((unsigned long)pow(P,i)*a)%(Q-1));
    if (b<a) return 0;
  }
  return 1;
}

// returns integer representation of the finite field element a
unsigned long integer_representation(ZZ_pE a) {
  ZZ_pX b;
  unsigned long c;
  conv(b,a); // converts a to a polynomial b
  unsigned long result;
  result = 0;
  for (int i=0; i<=deg(b); i++) {
    conv(c,coeff(b,i)); // converts the coefficient of X^i of b to an integer c
    result = result + (c*(unsigned long)pow(P,i));
  }
  return result;
}

// returns 1 if X^a is of DO-type. Otherwise, return 0
int is_quadratic(unsigned long a) {
  int vector_repr[N];
  int n_zero, n_one, n_two;

  // first convert a to vector
  unsigned long tmp_a=a;
  for (int i=0; i<N; i++) {
    vector_repr[i]=(tmp_a%P);
    tmp_a = (tmp_a/P);
  }

  n_zero = 0;
  n_one = 0;
  n_two = 0;
  for (int i=0; i<N; i++) {
    if (vector_repr[i]==0) {
      n_zero++;
    }
    if (vector_repr[i]==1) {
      n_one++;
    }
    if (vector_repr[i]==2) {
      n_two++;
    }
  }
  if ((n_zero==N-2) && (n_one==2)) return 1;
  if ((n_zero==N-1) && (n_two==1)) return 1;
  return 0;
}

// prints out all known planar monomials in this field
void print_known() {
  printf("Known planar monomials in this field:\n");
  for (int i=0; i<N; i++) {
    if (((P*N/(__gcd(i,(int)N)))%2)==1) {
      if (is_smallest_in_equivalence_class((unsigned long)pow(P,i)+1)) {
        printf("X^%lu\n", (unsigned long)pow(P,i)+1);
      }
    }
  }
  if (P==3) {
    for (int i=3; i<N; i++) {
      if (__gcd(i,2*(int)N)==1) {
        if (is_smallest_in_equivalence_class((unsigned long)((pow(P,i)+1)/2))) {
          printf("X^%lu\n", (unsigned long)((pow(P,i)+1)/2));
        }
      }
    }
  }
  printf("\n");
}

// evaluation of (x+1)^e - x^e
ZZ_pE eval(ZZ_pE x, unsigned long e) {
  ZZ_pE y,z;
  y = x+1;
  power(y,y,e);
  power(z,x,e);
  return (y-z);
}

// Return 1 if X^e is a candidate for being planar (in the sense that there is no contradiction after the collision search using N_SAMPLES samples)
// If 0 is returned, the monomial X^e cannot be planar
int is_planar_candidate(unsigned long e) {
  ZZ_pE x;
  ZZ_pE y;
  unsigned long i_x;
  unsigned long i_y;

  std::set<unsigned long> elements;
  std::set<unsigned long> queries;

  for (unsigned long i=0; i<N_SAMPLES; i++) {
    random(x); // choose random element x in GF(P^N)
    i_x = integer_representation(x);
    if (queries.find(i_x)==queries.end()) {
      y = eval(x,e); // compute y = (x+1)^e - x^e
      i_y = integer_representation(y);
      if (elements.find(i_y)!=elements.end()) return 0; // collision found, return 0
      elements.insert(i_y);
      queries.insert(i_x);
    }
  }
  return 1; // no collision found, X^e is a candidate for being planar
}

// returns 1 if x^k is planar on subfields p,p^2,p^3,p^4 and if conditions stated in the Coulter-Matthews paper are met. Returns 0 otherwise.
int is_planar_on_small_subfield(unsigned long k) {
  if (!(P==3)) {
    if (!((__gcd(k,Q-1)==2) && ((k%(P-1))==2))) return 0; // if those conditions are not met, exponent cannot be planar (see Coulter-Matthews paper)
  }
  if (P==3) {
    if (!((__gcd(k,Q-1)==2) && ((k%(2))==0))) return 0;
  }

  // we know that planar functions must be planar on every subfield. If F_Q contains a subfield F_P^2, F_P^3, or F_P^4, we use the characterization of planar exponents to check whether X^e can be planar
  // first, check conditions if F_P^2 is a subfield.
  if ((N%2)==0) {
      if (!(((k%(P2-1))==(2*P)) || ((k%(P2-1))==2))) return 0;
  }

  // then, check conditions if F_P^3 is a subfield
  if ((N%3)==0) {
      if (!( ((k%(P3-1))==(1+1)) || ((k%(P3-1))==(1+P)) || ((k%(P3-1))==(1+P2)) || ((k%(P3-1))==(P+P)) || ((k%(P3-1))==(P+P2)) || ((k%(P3-1))==(P2+P2)))) return 0;
  }

  // then, check conditions if F_P^4 is a subfield
  if (!(P==3)) {
    if ((N%4)==0) {
        if (!(((k%(P4-1))==(2*P)) || ((k%(P4-1))==2) || ((k%(P4-1))==(2*P2))  || ((k%(P4-1))==(2*P3)))) return 0;
      }
  }
  if ((P==3)) { // add the Coulter-Matthews cases to the previous cases
    if ((N%4)==0) {
        if (!(  ((k%(P4-1))==(P2+P+(2*1))) || ((k%(P4-1))==(P3+P2+(2*P))) || ((k%(P4-1))==(1+P3+(2*P2))) || ((k%(P4-1))==(P+1+(2*P3))) ||
                ((k%(P4-1))==(2*P)) || ((k%(P4-1))==2) || ((k%(P4-1))==(2*P2))  || ((k%(P4-1))==(2*P3)))) return 0;
      }
  }

  #ifdef USE_DIVISORS_5_TO_8
  if (N>=9) {
  if (!(P==3)) {
    if ((N%5)==0) {
      if (!( ((k%(P5-1))==((1+1)*1)) || ((k%(P5-1))==((1+P)*1)) || ((k%(P5-1))==((1+P2)*1)) || ((k%(P5-1))==((1+P3)*1)) || ((k%(P5-1))==((1+P4)*1)) ||
             ((k%(P5-1))==((1+1)*P)) || ((k%(P5-1))==((1+P)*P)) || ((k%(P5-1))==((1+P2)*P)) || ((k%(P5-1))==((1+P3)*P)) ||
             ((k%(P5-1))==((1+1)*P2)) || ((k%(P5-1))==((1+P)*P2)) || ((k%(P5-1))==((1+P2)*P2)) ||
             ((k%(P5-1))==((1+1)*P3)) || ((k%(P5-1))==((1+P)*P3)) ||
             ((k%(P5-1))==((1+1)*P4)) )) return 0;
    }
    if ((N%6)==0) {
      if (!( ((k%(P6-1))==((1+1)*1)) || ((k%(P6-1))==((1+P2)*1)) || ((k%(P6-1))==((1+P4)*1)) ||
             ((k%(P6-1))==((1+1)*P)) || ((k%(P6-1))==((1+P2)*P)) || ((k%(P6-1))==((1+P4)*P)) ||
             ((k%(P6-1))==((1+1)*P2)) || ((k%(P6-1))==((1+P2)*P2)) ||
             ((k%(P6-1))==((1+1)*P3)) || ((k%(P6-1))==((1+P2)*P3)) ||
             ((k%(P6-1))==((1+1)*P4)) ||
             ((k%(P6-1))==((1+1)*P5)) )) return 0;
    }
    if ((N%7)==0) {
      if (!( ((k%(P7-1))==((1+1)*1)) || ((k%(P7-1))==((1+P)*1)) || ((k%(P7-1))==((1+P2)*1)) || ((k%(P7-1))==((1+P3)*1)) || ((k%(P7-1))==((1+P4)*1)) || ((k%(P7-1))==((1+P5)*1)) || ((k%(P7-1))==((1+P6)*1)) ||
             ((k%(P7-1))==((1+1)*P)) || ((k%(P7-1))==((1+P)*P)) || ((k%(P7-1))==((1+P2)*P)) || ((k%(P7-1))==((1+P3)*P)) || ((k%(P7-1))==((1+P4)*P)) || ((k%(P7-1))==((1+P5)*P)) ||
             ((k%(P7-1))==((1+1)*P2)) || ((k%(P7-1))==((1+P)*P2)) || ((k%(P7-1))==((1+P2)*P2)) || ((k%(P7-1))==((1+P3)*P2)) || ((k%(P7-1))==((1+P4)*P2)) ||
             ((k%(P7-1))==((1+1)*P3)) || ((k%(P7-1))==((1+P)*P3)) || ((k%(P7-1))==((1+P2)*P3)) || ((k%(P7-1))==((1+P3)*P3)) ||
             ((k%(P7-1))==((1+1)*P4)) || ((k%(P7-1))==((1+P)*P4)) || ((k%(P7-1))==((1+P2)*P4)) ||
             ((k%(P7-1))==((1+1)*P5)) || ((k%(P7-1))==((1+P)*P5)) ||
             ((k%(P7-1))==((1+1)*P6)) )) return 0;
    }
    if ((N%8)==0) {
      if (!( ((k%(P8-1))==((1+1)*1)) || ((k%(P8-1))==((1+1)*P)) || ((k%(P8-1))==((1+1)*P2)) || ((k%(P8-1))==((1+1)*P3)) || ((k%(P8-1))==((1+1)*P4)) || ((k%(P8-1))==((1+1)*P5)) || ((k%(P8-1))==((1+1)*P6)) || ((k%(P8-1))==((1+1)*P7)) )) return 0;
    }
  }

  if ((P==3)) { // add the Coulter-Matthews cases to the previous cases
    if ((N%5)==0) {
      if (!( ((k%(P5-1))==(P2+P+(2*1))) || ((k%(P5-1))==(P3+P2+(2*P))) || ((k%(P5-1))==(P4+P3+(2*P2))) || ((k%(P5-1))==(1+P4+(2*P3))) || ((k%(P5-1))==(P+1+(2*P4))) ||
             ((k%(P5-1))==((1+1)*1)) || ((k%(P5-1))==((1+P)*1)) || ((k%(P5-1))==((1+P2)*1)) || ((k%(P5-1))==((1+P3)*1)) || ((k%(P5-1))==((1+P4)*1)) ||
             ((k%(P5-1))==((1+1)*P)) || ((k%(P5-1))==((1+P)*P)) || ((k%(P5-1))==((1+P2)*P)) || ((k%(P5-1))==((1+P3)*P)) ||
             ((k%(P5-1))==((1+1)*P2)) || ((k%(P5-1))==((1+P)*P2)) || ((k%(P5-1))==((1+P2)*P2)) ||
             ((k%(P5-1))==((1+1)*P3)) || ((k%(P5-1))==((1+P)*P3)) ||
             ((k%(P5-1))==((1+1)*P4)) )) return 0;
    }
    if ((N%6)==0) {
      if (!( ((k%(P6-1))==(P4+P3+P2+P+(2*1))) || ((k%(P6-1))==(P5+P4+P3+P2+(2*P))) || ((k%(P6-1))==(1+P5+P4+P3+(2*P2))) || ((k%(P6-1))==(P+1+P5+P4+(2*P3))) || ((k%(P6-1))==(P2+P+1+P5+(2*P4))) || ((k%(P6-1))==(P3+P2+P+1+(2*P5))) ||
             ((k%(P6-1))==((1+1)*1)) || ((k%(P6-1))==((1+P2)*1)) || ((k%(P6-1))==((1+P4)*1)) ||
             ((k%(P6-1))==((1+1)*P)) || ((k%(P6-1))==((1+P2)*P)) || ((k%(P6-1))==((1+P4)*P)) ||
             ((k%(P6-1))==((1+1)*P2)) || ((k%(P6-1))==((1+P2)*P2)) ||
             ((k%(P6-1))==((1+1)*P3)) || ((k%(P6-1))==((1+P2)*P3)) ||
             ((k%(P6-1))==((1+1)*P4)) ||
             ((k%(P6-1))==((1+1)*P5)) )) return 0;
    }
  }
  }
  #endif

  return 1;
}

int main(int argc, char* argv[])
{
  if (argc != 6) {
        	printf("Usage: %s <outfile> <P> <N> <start_exponent> <end_exponent>...\n", argv[0]);
        	return -1;
    	}

  FILE *fp = fopen(argv[1], "w");
  if (fp == NULL)
  {
      printf("Error opening file!\n");
  }

  P = (unsigned long)atoi(argv[2]);
  N = (unsigned long)atoi(argv[3]);
  Q = (unsigned long)pow(P,N);
  P2 = (unsigned long)pow(P,2);
  P3 = (unsigned long)pow(P,3);
  P4 = (unsigned long)pow(P,4);
  P5 = (unsigned long)pow(P,5);
  P6 = (unsigned long)pow(P,6);
  P7 = (unsigned long)pow(P,7);
  P8 = (unsigned long)pow(P,8);

  unsigned long start_exponent = (unsigned long)atoll(argv[4]);
  unsigned long end_exponent = (unsigned long)atoll(argv[5]);

  runtime = time(NULL);

  N_SAMPLES = 20*(unsigned long)sqrt(pow(P,N)); //this is for the collision search. By the birthday bound, we need roughly square root of P^N tests. We multiply by 20 to avoid false positives
  printf("Searching for planar monomials in the field of size %lu. Checking %lu samples for each exponent. Checking exponents in the range[%lu,%lu].\n\n",Q,N_SAMPLES,start_exponent,end_exponent);
  print_known();

   ZZ_p::init(ZZ(P)); // define GF(P)
   ZZ_pX P_IRR;
   BuildIrred(P_IRR, N); // generate the modulus as an irreducible polynomial of degree N over GF(P)
   ZZ_pE::init(P_IRR); // define GF(P^N)

   runtime = time(NULL);
   //unsigned long exponents_checked = 0;
   for (unsigned long e=start_exponent; e<end_exponent; e++) {
     if (is_smallest_in_equivalence_class(e) && is_planar_on_small_subfield(e)) {
       if (is_planar_candidate(e)) {
         printf("Candidate for planar monomial found: X^%lu, Quadratic: %d\n",e,is_quadratic(e));
         fprintf(fp,"Candidate for planar monomial found: X^%lu, Quadratic: %d\n",e,is_quadratic(e));
         fflush(stdout);
      }
    }
   }
   printf("Finished. Total running time: %li sec\n",time(NULL)-runtime);
   fclose(fp);
}
