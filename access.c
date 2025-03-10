
#include "hybridBV.h"
#include <time.h>
#include <sys/times.h>
#include <unistd.h>

uint64_t rnd (uint64_t m)

   { uint64_t r = rand();
     return r % m;
   }

uint64_t check (hybridBV B)

   { uint64_t l,r,s;
     if (B->type == tDynamic) {
        l = check(B->bv.dyn->left);
        r = check(B->bv.dyn->right);
	if (B->bv.dyn->leaves != l+r) {
	   printf("fallo, %li != %li + %li\n",B->bv.dyn->leaves,l,r);
	   exit(1);
	   }
        return l+r;
	}
     return hybridLeaves(B);
   }

void main (int argc, char **argv)

   { hybridBV B;
     uint64_t n,m,i,o,u;
     uint64_t *data;
     struct tms t1,t2;
     float alphaUpd;
     int breve;

     if (argc < 4)
        { fprintf(stderr,"Usage: %s <log_2 n> <alpha> <1/q> [<factor>]\n"
          "Creates a bitvector of length n and applies alpha*n ops on it,\n"
          "where a fraction 1/q of them are indels (50/50 in probability)\n"
          "and the others are accesses, all at random positions.\n"
          "Nodes containing t bits are flattened after receiving factor*t queries.\n",
                  argv[0]);
          exit(1);
        }

     srand(time(NULL)); 

     u = n = (((uint64_t)1) << atoi(argv[1]));
     m = n * atoi(argv[2]);
     alphaUpd = atof(argv[3]);
     if (argc > 4) Theta = atof(argv[4]);
     breve = (argc > 5) && !strcmp(argv[5],"-");

     data = (uint64_t*)myalloc(n/8);
     for (i=0;i<n/w;i++)
         data[i] = rand() | (((uint64_t)rand()) << 32);

     times(&t1);

     B = hybridCreateFrom(data,n);

     for (i=0;i<m;i++)
	{ if (rnd(1000000000) < alphaUpd * 1000000000)
             { if (rnd(2)==0) hybridInsert(B,rnd(++n),rnd(2));
               else hybridDelete(B,rnd(n--));
             } 
          else  
             { if (rnd(1)==0) hybridAccess(B,rnd(n));
               else hybridSelect(B,1+rnd(hybridOnes(B)));
	     }
        }
     times(&t2);
     if (breve) printf ("%f + ",
	     (t2.tms_utime-t1.tms_utime)/(float)sysconf(_SC_CLK_TCK)*1000000/(float)m);
     else printf("Time per operation in microseconds: %f\n",
	     (t2.tms_utime-t1.tms_utime)/(float)sysconf(_SC_CLK_TCK)*1000000/(float)m);
  /*   printf ("n = %li, m = %li, alphaUpd = %f, utime = %f microseconds\n",
	     u,m,alphaUpd,
	     (t2.tms_utime-t1.tms_utime)/(float)sysconf(_SC_CLK_TCK)*1000000/(float)m);
*/
   //  printf("   %li hojas\n",check (B));

     hybridDestroy(B);
     exit(0);
   }
