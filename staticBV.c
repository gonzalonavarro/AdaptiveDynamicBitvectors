
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

#include "staticBV.h"

#define K 4 // block length is w*K

#define w16 (8*sizeof(uint16_t)) // superblock length is 2^w16

	// preprocesses for rank, with parameter K

static void staticPreprocess (staticBV B)

    { uint64_t i,n;
      uint64_t sacc,acc;
      uint64_t last,word;
      n = B->size;
      if (n == 0) return;
      B->B = (uint16_t*)myalloc(((n+K*w-1)/(K*w))*sizeof(uint16_t));
      B->S = (uint64_t*)myalloc(((n+(1<<w16)-1)/(1<<w16))*sizeof(uint64_t));
      sacc = acc = 0;
      i = 0;
      while (i<(n+w-1)/w)
          { uint64_t top = min((n+w-1)/w,i+(1<<w16)/w);
            sacc += acc; acc = 0;
            B->S[(i*w) >> w16] = sacc;
            while (i<top)
                { if (i%K == 0) B->B[i/K] = acc;
                  acc += popcount(B->data[i]);
                  i++;
                }
          }
      B->ones = staticRank(B,n-1);
    } 

	// converts a bit array into a bitvector of n bits
        // data is pointed to and will be freed 

staticBV staticCreateFrom (uint64_t *data, uint64_t n)

    { staticBV B;
      B = (staticBV)myalloc(sizeof(struct s_staticBV));
      B->size = n;
      if (n == 0) B->data = NULL;
      else B->data = data;
      B->S = NULL;
      B->B = NULL;
      staticPreprocess(B);
      return B;
    }

	// destroys B, frees data

void staticDestroy (staticBV B)

    { if (B != NULL) 
         { myfree(B->data);
           myfree(B->S); myfree(B->B);
      	   myfree(B);
	 }
    }

        // writes B's data to file, which must be opened for writing 

void staticSave (staticBV B, FILE *file)

   { if (B->size != 0)
        fwrite (B->data,sizeof(uint64_t),(B->size+w-1)/w,file);
   }

        // loads staticBV's data from file, which must be opened for reading
	// size is the number of bits

staticBV staticLoad (FILE *file, uint64_t size)

    { staticBV B;
      B = (staticBV)myalloc(sizeof(struct s_staticBV));
      B->size = size;
      if (size == 0) B->data = NULL;
      else { B->data = (uint64_t*)myalloc(((B->size+w-1)/w)*sizeof(uint64_t));
	     fread (B->data,sizeof(uint64_t),(B->size+w-1)/w,file);
	   }
      B->S = NULL;
      B->B = NULL;
      staticPreprocess(B);
      return B;
    }

	// data of staticBV

extern inline uint64_t *staticBits (staticBV B)

    { return B->data;
    }

	// staticBV size in w-bit words

uint64_t staticSpace (staticBV B)

    { uint64_t space = sizeof(struct s_staticBV)*8/w;
      if (B == NULL) return 0;
      if (B->data != NULL) space += (B->size+w-1)/w;
      if (B->B != NULL) space += ((B->size+K*w-1)/(K*w))/(w/w16);
      if (B->S != NULL) space += (B->size+(1<<w16)-1)/(1<<w16);
      return space;
    }

        // gives bit length

extern inline uint64_t staticLength (staticBV B)

    { return B->size;
    }

        // gives number of ones

extern inline uint64_t staticOnes (staticBV B)

    { return B->ones;
    }

	// access B[i], assumes i is right

extern inline uint staticAccess (staticBV B, uint64_t i)

    { return (B->data[i/w] >> (i%w)) & 1;
    }

        // read bits [i..i+l-1] onto D[j..], assumes it is right
        
extern inline void staticRead (staticBV B, uint64_t i, uint64_t l, 
			       uint64_t * D, uint64_t j)

    { copyBits(D,j,B->data,i,l);
    }

	// computes rank(B,i), zero-based, assumes i is right

extern inline uint64_t staticRank (staticBV B, uint64_t i)

    { uint64_t b,sb;
      uint64_t rank;
      sb = i/(K*w);
      rank = B->S[i>>w16] + B->B[sb];
      sb *= K;
      for (b=sb;b<i/w;b++) rank += popcount(B->data[b]);
      return rank + popcount(B->data[b] & (((uint64_t)~0) >> (w-1-(i%w))));
    }

	// computes rank(B,i), zero-based, assumes i is right

extern inline uint64_t staticRank0 (staticBV B, uint64_t i)

    { return i + 1 - staticRank(B,i);
    }

        // computes select_1(B,j), zero-based, assumes j is right

extern uint64_t staticSelect (staticBV B, uint64_t j)

    { int64_t i,d,b;
      uint p;
      uint64_t word,s,m,n;
      n = B->size;
      s = (n+(1<<w16)-1)/(1<<w16);
	// interpolation: guess + exponential search
      i = ((uint64_t)(j * (n / (float)B->ones))) >> w16;
      if (i == s) i--;
      if (B->S[i] < j)
	 { d = 1;
	   while ((i+d < s) && (B->S[i+d] < j))
	      { i += d; d <<= 1; }
	   d = min(s,i+d); // now d is the top of the range
	   while (i+1<d)
	      { m = (i+d)>>1;
	        if (B->S[m] < j) i = m; else d = m;
	      }
	 }
      else
	 { d = 1;
	   while ((i-d >= 0) && (B->S[i-d] >= j))
	      { i -= d; d <<= 1; }
	   d = max(0,i-d); // now d is the bottom of the range
	   while (d+1<i)
	      { m = (i+d)>>1;
	        if (B->S[m] < j) d = m; else i = m;
	      }
	   i--;
	 }
	// now the same inside the superblock
      j -= B->S[i]; // what remains to be found inside the superblock
      p = i < s-1 ? B->S[i+1]-B->S[i] : B->ones-B->S[i];
      b = (i << w16) / (w*K);
      s = min(b+(1<<w16)/(w*K),(n+w*K-1)/(w*K));
      i = b + ((j * (s-b)*(w*K) / (float)p)) / (w*K);
      if (i == s) i--;
      if (B->B[i] < j)
	 { d = 1;
	   while ((i+d < s) && (B->B[i+d] < j))
	      { i += d; d <<= 1; }
	   d = min(s,i+d); // now d is the top of the range
	   while (i+1<d)
	      { m = (i+d)>>1;
	        if (B->B[m] < j) i = m; else d = m;
	      }
	 }
      else
	 { d = 1;
	   while ((i-d >= b) && (B->B[i-d] >= j))
	      { i -= d; d <<= 1; }
	   d = max(b,i-d); // now d is the bottom of the range
	   while (d+1<i)
	      { m = (i+d)>>1;
	        if (B->B[m] < j) d = m; else i = m;
	      }
	   i--;
	 }
	// now it's confined to K blocks
      j -= B->B[i];
      i *= K;
      while ((i+1)*w < n)
	{ p = popcount(B->data[i]);
	  if (p >= j) break;
	  j -= p;
	  i++;
	}
      word = B->data[i];
      i *= w;
/* this was actually slower
     uint len = (8*sizeof(word)) >> 1;
     while (len)
        { p = popcount(word & ((((uint64_t)1) << len)-1));
          if (p < j)
             { word >>= len;
               j -= p;
               i += len;
             }
          len >>= 1;
        }
     return i;
*/
      while (1)
	{ j -= word & 1;
	  if (j == 0) return i;
	  word >>= 1;
	  i++;
	}
    }

        // computes select_0(B,j), zero-based, assumes j is right

extern uint64_t staticSelect0 (staticBV B, uint64_t j)

    { int64_t i,d,b;
      uint p;
      uint64_t word,s,m,n,sb;
      n = B->size;
      s = (n+(1<<w16)-1)/(1<<w16);
	// interpolation: guess + exponential search
      i = ((uint64_t)(j * (n / (float)(n-B->ones)))) >> w16;
      if (i == s) i--;
      if (i*(1 << w16) - B->S[i] < j)
	 { d = 1;
	   while ((i+d < s) && ((i+d)*(1 << w16) - B->S[i+d] < j))
	      { i += d; d <<= 1; }
	   d = min(s,i+d); // now d is the top of the range
	   while (i+1<d)
	      { m = (i+d)>>1;
	        if (m*(1 << w16) - B->S[m] < j) i = m; else d = m;
	      }
	 }
      else
	 { d = 1;
	   while ((i-d >= 0) && ((i-d)*(1 << w16) - B->S[i-d] >= j))
	      { i -= d; d <<= 1; }
	   d = max(0,i-d); // now d is the bottom of the range
	   while (d+1<i)
	      { m = (i+d)>>1;
	        if (m*(1 << w16) - B->S[m] < j) d = m; else i = m;
	      }
	   i--;
	 }
	// now the same inside the superblock
      j -= i*(1 << w16) - B->S[i]; // what remains to be found inside superblock
      p = i < s-1 ? (1<<w16)-(B->S[i+1]-B->S[i]) : 
		    (B->size-i*(1<<w16))-(B->ones-B->S[i]);
      b = (i << w16) / (w*K);
      s = min(b+(1<<w16)/(w*K),(n+w*K-1)/(w*K));
      i = b + ((j * (s-b)*(w*K) / (float)p)) / (w*K);
      if (i == s) i--;
      if ((i-b)*w*K - B->B[i] < j)
	 { d = 1;
	   while ((i+d < s) && (((i-b)+d)*w*K - B->B[i+d] < j))
	      { i += d; d <<= 1; }
	   d = min(s,i+d); // now d is the top of the range
	   while (i+1<d)
	      { m = (i+d)>>1;
	        if ((m-b)*w*K - B->B[m] < j) i = m; else d = m;
	      }
	 }
      else
	 { d = 1;
	   while ((i-d >= b) && (((i-b)-d)*w*K - B->B[i-d] >= j))
	      { i -= d; d <<= 1; }
	   d = max(b,i-d); // now d is the bottom of the range
	   while (d+1<i)
	      { m = (i+d)>>1;
	        if ((m-b)*w*K - B->B[m] < j) d = m; else i = m;
	      }
	   i--;
	 }
	// now it's confined to K blocks
      j -= (i-b)*w*K - B->B[i];
      i *= K;
      while ((i+1)*w < n)
	{ p = popcount(~B->data[i]);
	  if (p >= j) break;
	  j -= p;
	  i++;
	}
      word = ~B->data[i];
      i *= w;
      while (1)
	{ j -= word & 1;
	  if (j == 0) return i;
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

int64_t staticNext (staticBV B, uint64_t i)

    { uint64_t p,b,sb;
      uint64_t word,rank;

      p = i/w;
      word = B->data[p] & ((~(uint64_t)0)<<(i%w));
      if ((p+1) * w > B->size) 
         word &= (((uint64_t)1) << (B->size % w)) - 1;
      if (word) // a likely case, solve faster
         return p*w + decode[(0x03f79d71b4ca8b09 * (word & -word))>>58];
	// search within block
      b = min((p/K+2)*K,1+(B->size-1)/w); // scan at least 2 blocks (a full one)
      p++;
      while (p < b)
         { word = B->data[p++];
	   if (word) break;
	 }
      if (word)
         { if (p*w > B->size)
              word &= (((uint64_t)1) << (B->size % w)) - 1;
           if (!word) return -1; // end of bitvector
           return (p-1)*w + decode[(0x03f79d71b4ca8b09 * (word & -word))>>58];
	 }
      else if (p == 1+(B->size-1)/w) return -1; // end of bitvector
	// reduce to select
      sb = p/K-1;
      rank = B->S[(sb*K*w)>>w16] + B->B[sb];
      if (rank == B->ones) return -1;
      return staticSelect(B,rank+1);
    }
		
        // computes next_0(B,i), zero-based and including i
        // returns -1 if no answer

int64_t staticNext0 (staticBV B, uint64_t i)

    { uint64_t p,b,sb;
      uint64_t word,rank;

      p = i/w;
      word = ~B->data[p] & ((~(uint64_t)0)<<(i%w));
      if (word) // a likely case, solve faster
         return p*w + decode[(0x03f79d71b4ca8b09 * (word & -word))>>58];
	// search within block
      b = min((p/K+2)*K,1+(B->size-1)/w); // scan at least 2 blocks (a full one)
      p++;
      while (p < b)
         { word = ~B->data[p++];
	   if (word) break;
	 }
      if (word)
         { if (p*w > B->size)
              word &= (((uint64_t)1) << (B->size % w)) - 1;
           if (!word) return -1; // end of bitvector
           return (p-1)*w + decode[(0x03f79d71b4ca8b09 * (word & -word))>>58];
	 }
      else if (p == 1+(B->size-1)/w) return -1; // end of bitvector
	// reduce to select
      sb = p/K-1;
      rank = sb*K*w - (B->S[(sb*K*w)>>w16] + B->B[sb]);
      if (rank == B->size - B->ones) return -1;
      return staticSelect0(B,rank+1);
    }
		
