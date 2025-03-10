
/*

HybridBV -- an implementation of adaptive dynamic bitvectors. 
Copyright (C) 2024-current_year Gonzalo Navarro

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

Author's contact: Gonzalo Navarro, Dept. of Computer Science, University of
Chile. Beauchef 851, Santiago, Chile. gnavarro@dcc.uchile.cl

*/

	// supports leaf bitvectors of size up to 2^32-1

#include "leafBV.h"

const int MaxBlockWords = 128; // b value in words: maximum leaf size
const float Gamma = 0.75; // new blocks try to be this fraction full

       // size that a newly created leaf should have

extern inline uint leafNewSize(void)

   { return MaxBlockWords * Gamma;
   }

        // max leaf size

extern inline uint leafMaxSize(void)

   { return MaxBlockWords;
   }

	// creates an empty leafBV 

leafBV leafCreate (void)

   { leafBV B = (leafBV)myalloc(sizeof(struct s_leafBV));
     B->size = 0;
     B->ones = 0;
     B->data = (uint64_t*)myalloc(MaxBlockWords * sizeof(uint64_t));
     return B;
   }

	// converts a bit array into a leafBV of n bits
	// data is copied but freed only if freeit

leafBV leafCreateFrom (uint64_t *data, uint n, int freeit)

   { uint i,nb;
     leafBV B;
     if (n == 0) return leafCreate();
     B = (leafBV)myalloc(sizeof(struct s_leafBV));
     B->size = n;
     B->data = (uint64_t*)myalloc(MaxBlockWords * sizeof(uint64_t));
     nb = (n+7)/8;
     memcpy(B->data,data,nb);
     nb = (n+w-1)/w;
     if (freeit) free(data);
     if (n % w) B->data[nb-1] &= (((uint64_t)1) << (n % w)) - 1;
     B->ones = 0;
     for (i=0;i<nb;i++) B->ones += popcount(B->data[i]);
     return B;
   }

	// destroys B, frees data 

void leafDestroy (leafBV B)

   { myfree(B->data);
     myfree(B);
   }

	// saves leaf data to file, which must be opened for writing

void leafSave (leafBV B, FILE *file)

   { if (B->size != 0)
        myfwrite (B->data,sizeof(uint64_t),(B->size+w-1)/w,file);
   }

        // loads leaf data from file, which must be opened for reading
        // size is the number of bits

leafBV leafLoad (FILE *file, uint size)

    { uint64_t *data = (uint64_t*)myalloc(((size+w-1)/w)*sizeof(uint64_t));
      myfread (data,sizeof(uint64_t),(size+w-1)/w,file);
      return leafCreateFrom(data,size,1);
    }

	// gives (allocated) space of B in w-bit words

uint leafSpace (leafBV B)

   { return (sizeof(struct s_leafBV)+sizeof(uint64_t)-1)/sizeof(uint64_t) 
	    + MaxBlockWords;
   }

	// gives bit length

extern inline uint leafLength (leafBV B)

   { return B->size;
   }

	// gives number of 1s

extern inline uint leafOnes (leafBV B)

   { return B->ones;
   }

       // sets value for B[i]= (v != 0), assumes i is right
        // returns difference in 1s

int leafWrite (leafBV B, uint i, uint v)

   { uint64_t one = ((uint64_t)1) << (i%w);
     if (v) {
	if (!(B->data[i/w] & one)) {
	   B->data[i/w] |= one;
	   B->ones++;
	   return 1;
	   }
	}
     else {
	if (B->data[i/w] & one) {
	   B->data[i/w] &= ~one;
	   B->ones--;
	   return -1;
	   }
	}
     return 0;
   }

        // inserts v at B[i], assumes i is right and that insertion is possible

void leafInsert (leafBV B, uint i, uint v)

   { uint nb = ++B->size/w;
     uint ib = i/w;
     int b;

     for (b=nb;b>ib;b--)
	 B->data[b] = (B->data[b] << 1) | (B->data[b-1] >> (w-1));
     if ((i+1)%w)
          B->data[ib] = (B->data[ib] & ((((uint64_t)1) << (i%w)) - 1)) |
		        (((uint64_t)v) << (i%w)) |
		        ((B->data[ib] << 1) & (~((uint64_t)0) << ((i+1)%w)));
     else B->data[ib] = (B->data[ib] & ((((uint64_t)1) << (i%w)) - 1)) |
                        (((uint64_t)v) << (i%w));
     B->ones += v;
   }

        // deletes B[i], assumes i is right
        // returns difference in 1s

int leafDelete (leafBV B, uint i)

   { uint nb = B->size--/w;
     uint ib = i/w;
     int b;
     int v = (B->data[ib] >> (i%w)) & 1;

     B->data[ib] = (B->data[ib] & ((((uint64_t)1) << (i%w)) - 1)) |
		   ((B->data[ib] >> 1) & (~((uint64_t)0) << (i%w)));
     for (b=ib+1;b<=nb;b++) {
	 B->data[b-1] |= B->data[b] << (w-1);
	 B->data[b] >>= 1;
	 }
     B->ones -= v;
     return -v;
   }

	// access B[i], assumes i is right

extern inline uint leafAccess (leafBV B, uint i)

   { return (B->data[i/w] >> (i%w)) & 1;
   }

        // read bits [i..i+l-1], onto D[j...]

void leafRead (leafBV B, uint i, uint l, uint64_t *D, uint64_t j)

   { copyBits(D,j,B->data,i,l);
   }

	// computes rank(B,i), zero-based, assumes i is right

uint leafRank (leafBV B, uint i)

   { int p,ib;
     uint ones = 0;
     ib = ++i/w;
     for (p=0;p<ib;p++) ones += popcount(B->data[p]);
     if (i%w) ones += popcount(B->data[p] & ((((uint64_t)1)<<(i%w))-1));
     return ones;
   }

        // computes select_1(B,j), zero-based, assumes j is right

uint leafSelect (leafBV B, uint j)

   { uint p,i,pc;
     uint64_t word;
     uint ones = 0;
     p = 0;
     while (1)
	{ word = B->data[p];
	  pc = popcount(word);
	  if (ones+pc >= j) break;
	  ones += pc; 
	  p++;
	}
     i = p*w;
/* this was actually slower
     uint len = (8*sizeof(word)) >> 1;
     while (len)
	{ pc = popcount(word & ((((uint64_t)1) << len)-1)); 
          if (ones + pc < j)
	     { word >>= len;
	       ones += pc;
	       i += len;
	     }
	  len >>= 1;
	}
     return i;
*/
     i = p*w;
     while (1)
	{ ones += word & 1; 
	  if (ones == j) return i;
	  word >>= 1; 
	  i++; 
	}
   }

        // computes select_1(B,j), zero-based, assumes j is right

uint leafSelect0 (leafBV B, uint j)

   { uint p,i,pc;
     uint64_t word;
     uint ones = 0;
     p = 0;
     while (1)
	{ word = ~B->data[p];
	  pc = popcount(word);
	  if (ones+pc >= j) break;
	  ones += pc; 
	  p++;
	}
     i = p*w;
     while (1)
	{ ones += word & 1; 
	  if (ones == j) return i;
	  word >>= 1; 
	  i++; 
	}
   }

        // computes next_1(B,i), zero-based and including i
	// returns -1 if no answer

	// trick for lowest 1 in a 64-bit word
static int decode[64] = {
       0, 1,56, 2,57,49,28, 3,61,58,42,50,38,29,17, 4,
      62,47,59,36,45,43,51,22,53,39,33,30,24,18,12, 5,
      63,55,48,27,60,41,37,16,46,35,44,21,52,32,23,11,
      54,26,40,15,34,20,31,10,25,14,19, 9,13, 8, 7, 6 };

int leafNext (leafBV B, uint i)

   { uint p;
     uint64_t word;
     p = i/w;
     word = B->data[p] & ((~(uint64_t)0)<<(i%w));
     while ((++p*w <= B->size) && !word)
	word = B->data[p];
     if (p*w > B->size)
	word &= (((uint64_t)1) << (B->size % w)) - 1;
     if (!word) return -1;
     return (p-1)*w + decode[(0x03f79d71b4ca8b09 * (word & -word))>>58];
   }

        // computes next_0(B,i), zero-based and including i
	// returns -1 if no answer

int leafNext0 (leafBV B, uint i)

   { uint p;
     uint64_t word;
     p = i/w;
     word = (~B->data[p]) & ((~(uint64_t)0)<<(i%w));
     while ((++p*w <= B->size) && word)
	word = ~B->data[p];
     if (p*w > B->size)
	word &= (((uint64_t)1) << (B->size % w)) - 1;
     if (!word) return -1;
     return (p-1)*w + decode[(0x03f79d71b4ca8b09 * (word & -word))>>58];
   }
