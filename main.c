
#include "hybridBV.h"
#include "hybridId.h"
#include <time.h>
#include <sys/times.h>
#include <unistd.h>

// #define ADV
// #define ADV2
// #define BASIC
// #define BASICID
// #define ADVID
// #define ADVID2
// #define WORSTCASE
#define NEXT

uint64_t rnd (uint64_t m)

   { uint64_t r = rand();
     return r % m;
   }

#define alphaUpd 1.0

void main (void)

   { hybridBV B;
     hybridId I;
     uint64_t n,m,i,o,u;
     uint64_t *data;
     staticBV S;
     struct tms t1,t2;
     int64_t j,k;

     srand(time(NULL)); 

#ifdef NEXT

     n = 1024*1024*w;
     m = n * 0.1;

     data = (uint64_t*)calloc(n/8,1);
     for (i=0;i<m;i++)
         { j = rand() % n;
	   data[j/w] |= ((uint64_t)1) << (j%w);
	 }
     // for 0: for (i=0;i<n/w;i++) data[i] = ~data[i];

     B = hybridCreateFrom(data,n);

     times(&t1);
     
     j = 0; k = 1;
     while (1)
        { j = hybridNext(B,j);
	  if (j == -1) break;
          if (hybridSelect(B,k++) != j) printf("Mal!\n");
	  j++;
	}
     if (hybridRank(B,n-1) != k-1) printf("Mal\n");

     times(&t2);
     printf ("utime = %.2f microseconds\n",(t2.tms_utime-t1.tms_utime)/(float)sysconf(_SC_CLK_TCK)*1000000/(float)(n/leafIdNewSize(25)));

     exit(1); 

#endif

#ifdef WORSTCASE

     u = n = 1024*1024*w;
     m = 1024*1024*10;

     data = (uint64_t*)myalloc(n*sizeof(uint64_t));
     for (i=0;i<n;i++)
         data[i] = i;

     times(&t1);

     I = hybridIdCreateFrom64(data,n,25);

     for (i=0;i<n;i++)
         hybridIdAccess(I,i);
        
     for (i=0;i<n;i+=leafIdNewSize(25))
	   hybridIdWrite(I,i,0);

     times(&t2);
     printf ("utime = %.2f microseconds\n",(t2.tms_utime-t1.tms_utime)/(float)sysconf(_SC_CLK_TCK)*1000000/(float)(n/leafIdNewSize(25)));

     hybridIdDestroy(I);

#endif

#ifdef ADVID

     u = n = 1024*1024*w;
     m = 1024*1024*10;

     data = (uint64_t*)myalloc(n*sizeof(uint64_t));
     for (i=0;i<n;i++)
         data[i] = i;

     times(&t1);

     I = hybridIdCreateFrom64(data,n,25);

     for (i=0;i<m;i++)
	{ if (rnd(1000000000) < alphaUpd * 1000000000)
             { if (rnd(2)==0) hybridIdInsert(I,rnd(++n),rnd(u));
               else hybridIdDelete(I,rnd(n--));
             } 
          else  
             { hybridIdAccess(I,rnd(n));
	     }
        }
     times(&t2);
     printf ("utime = %.2f microseconds\n",(t2.tms_utime-t1.tms_utime)/(float)sysconf(_SC_CLK_TCK)*1000000/(float)m);

     hybridIdDestroy(I);

#endif

#ifdef ADVID2

     u = n = 1024*1024*w;
     m = 1024*1024*10;

     data = (uint64_t*)myalloc(n*sizeof(uint64_t));
     for (i=0;i<n;i++)
         data[i] = i;

     times(&t1);

     I = hybridIdCreateFrom64(data,n,25);

     for (i=0;i<m;i++)
	{ if (rnd(1000000000) < alphaUpd * 1000000000)
             { if (rnd(2)==0) hybridIdInsert(I,0*++n,rnd(u));
               else hybridIdDelete(I,--n); 
             } 
          else  
             { hybridIdAccess(I,rnd(n));
	     }
        }
     times(&t2);
     printf ("utime = %.2f microseconds\n",(t2.tms_utime-t1.tms_utime)/(float)sysconf(_SC_CLK_TCK)*1000000/(float)m);

     hybridIdDestroy(I);

#endif

#ifdef ADV2

     n = 1024*1024*w;
     m = 1024*1024*10;

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
	     { o = 1+rnd(hybridOnes(B));
	       if (hybridRank0(B,hybridSelect0(B,o)) != o)
	           { printf ("Mal!\n"); }
	     }
        }
     times(&t2);
     printf ("utime = %.2f microseconds\n",(t2.tms_utime-t1.tms_utime)/(float)sysconf(_SC_CLK_TCK)*1000000/(float)m);

     hybridDestroy(B);

#endif

#ifdef ADV

     n = 1024*1024*w;
     m = 1024*1024*10;

     data = (uint64_t*)myalloc(n/8);
     for (i=0;i<n/w;i++)
         data[i] = rand() | (((uint64_t)rand()) << 32);

     times(&t1);

     B = hybridCreateFrom(data,n);

     for (i=0;i<m;i++)
	{ if (rnd(1000000000) < alphaUpd * 1000000000)
             { if (rnd(2)==0) hybridInsert(B,0*++n,rnd(2)); // hybridInsert(B,rnd(++n),rnd(2));
               else hybridDelete(B,--n); // hybridDelete(B,rnd(n--));
             } 
          else  
	     { // o = 1+rnd(hybridOnes(B));
	       // if (hybridRank(B,hybridSelect(B,o)) != o)
	       //    { printf ("Mal!\n"); }
	       // continue;
               if (rnd(2)==0) hybridRank(B,rnd(n));
	       else hybridSelect(B,1+rnd(hybridOnes(B)));
	     }
        }
     times(&t2);
     printf ("utime = %.2f microseconds\n",(t2.tms_utime-t1.tms_utime)/(float)sysconf(_SC_CLK_TCK)*1000000/(float)m);

     hybridDestroy(B);

#endif

#ifdef BASIC

     B = hybridCreate();

     printf("Inserting %li 10s\n",leafMaxSize()*w);
     for (i=0;i<leafMaxSize()*w;i++)
         { hybridInsert(B,0,0);
           hybridInsert(B,0,1);
	 }

     printf("Reading many times to make it static\n");
     for (i=0;i< leafMaxSize()*w*10;i++)
         hybridAccess(B,i);

     printf("B[0..%li] (%li ones) = ",hybridLength(B)-1,hybridOnes(B));
     for (i=0;i<hybridLength(B);i++)
	 printf("%i",hybridAccess(B,i));
     printf("\n");

     printf("Ranks at 10: ");
     for (i=0;i<hybridLength(B);i+=10)
	 printf("%li,",hybridRank(B,i));
     printf("\n");

     printf("Deleting all the 0s\n");
     for (i=0;i<leafMaxSize()*w;i++)
         { hybridDelete(B,i+1);
	 }
     printf("B[0..%li] (%li ones) = ",hybridLength(B)-1,hybridOnes(B));
     for (i=0;i<hybridLength(B);i++)
	 printf("%i",hybridAccess(B,i));
     printf("\n");

     printf("Writing 0s every other position\n");
     for (i=0;i<hybridLength(B);i+=2)
         { hybridWrite(B,i,0);
	 }
     printf("B[0..%li] (%li ones) = ",hybridLength(B)-1,hybridOnes(B));
     for (i=0;i<hybridLength(B);i++)
	 printf("%i",hybridAccess(B,i));
     printf("\n");

     printf("Deleting half of the bits\n");
     for (i=0;i<hybridLength(B)/2;i++)
         { hybridDelete(B,0);
	 }
     printf("B[0..%li] (%li ones) = ",hybridLength(B)-1,hybridOnes(B));
     for (i=0;i<hybridLength(B);i++)
	 printf("%i",hybridAccess(B,i));
     printf("\n");

     hybridDestroy(B);

#endif

#ifdef BASICID

     I = hybridIdCreate(10);

     printf("Inserting 1 to 1000\n");
     for (i=0;i<1000;i++)
         { hybridIdInsert(I,i,i+1);
	 }

     printf("Reading...\n");
     for (i=0;i<1000;i++)
         { printf("%li ",o = hybridIdAccess(I,i));
	   if (o != i+1) printf ("Mal\n");
	 }
     printf("\n");

     printf("Changing i to 1001-i\n");
     for (i=0;i<1000;i++)
         { hybridIdWrite(I,i,1000-i);
	 }

     printf("Reading...\n");
     for (i=0;i<1000;i++)
         { printf("%li ",o = hybridIdAccess(I,i));
	   if (o != 1000-i) printf ("Mal\n");
	 }
     printf("\n");

     printf("Deleting 101 to 900\n");
     for (i=100;i<900;i++)
         { hybridIdDelete(I,100);
	 }

     printf("Reading...\n");
     for (i=0;i<200;i++)
         { printf("%li ",o = hybridIdAccess(I,i));
	 }
     printf("\n");

     hybridIdDestroy(I);

#endif

   }
