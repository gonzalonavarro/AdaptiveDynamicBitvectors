
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

#include "basics.h"

void *myalloc (size_t n)

   { void *p;
     if (n == 0) return NULL;
     p = malloc(n);
     if (p == NULL) 
        { fprintf(stderr,"Error: malloc of %li bytes returned null\n",n);
          exit(1);
        }
     return p;
   }

void *mycalloc (size_t n, size_t s)

   { void *p;
     if (n == 0) return NULL;
     p = calloc(n,s);
     if (p == NULL) 
        { fprintf(stderr,"Error: calloc of %li bytes returned null\n",n);
          exit(1);
        }
     return p;
   }

void *myrealloc (void *p, size_t n)

   { if (p == NULL) return myalloc(n);
     if (n == 0) return NULL;
     p = realloc(p,n);
     if (p == NULL) 
        { fprintf(stderr,"Error: realloc of %li bytes returned null\n",n);
	  exit(1);
	}
     return p;
   }

void myfree (void *p)

  { if (p != NULL) free(p);
  }

void myfread(void *ptr, size_t size, size_t nmemb, FILE *stream)

  { if (fread(ptr,size,nmemb,stream) != nmemb)
       { fprintf(stderr,"Error: fread of %li bytes failed\n",nmemb*size);
	 exit(1);
       }
  }

void myfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)

  { if (fwrite(ptr,size,nmemb,stream) != nmemb)
       { fprintf(stderr,"Error: fwrite of %li bytes failed\n",nmemb*size);
	 exit(1);
       }
  }

uint numbits (uint n)

   { uint bits = 0;
     while (n)
        { n = n>>1; bits++; }
     return bits ? bits : 1;
   }

extern inline uint popcount (uint64_t y)

    { y -= ((y >> 1) & 0x5555555555555555ull);
      y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
      return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
    }

        // copies len bits starting at *src + psrc
        // to tgt from bit position ptgt
        // WARNING: writes some extra bits after target (but not more words)

void copyBits (uint64_t *tgt, uint64_t ptgt,
               uint64_t *src, uint64_t psrc, uint64_t len)

   { uint64_t old,mask;

        // easy cases if they are similarly misaligned
        // I didn't align to 8 due to possible endianness issues
     tgt += ptgt/w; ptgt %= w;
     src += psrc/w; psrc %= w;
     mask = (((uint64_t)1) << ptgt)-1;
     if (ptgt == psrc)
        { if (ptgt != 0)
             { *tgt = (*tgt & mask) + (*src & ~mask);
               *tgt++; *src++; len -= w-ptgt;
             }
          memcpy (tgt,src,((len+w-1)/w)*sizeof(uint64_t));
          return;
        }
        // general case, we first align the source
     if (ptgt < psrc) // first word from src fits in ptgt
        { *tgt = (*tgt & mask) + ((*src++ >> (psrc-ptgt)) & ~mask);
          ptgt += w-psrc;
        }
     else
        { *tgt = (*tgt & mask) + ((*src << (ptgt-psrc)) & ~mask);
          if (len <= w-ptgt) return; // before overflowing to the next word
          ptgt -= psrc;
          *++tgt = *src++ >> (w-ptgt); 
        }
     if (len <= w-psrc) return;
     len -= w-psrc;
        // now src is aligned, copy all the rest
     mask = (((uint64_t)1) << ptgt)-1;
     old = *tgt & mask;
     len += w; // cannot write len >= 0 as it is unsigned
     while (len > w)
        { *tgt++ = old + (*src << ptgt);
          old = *src++ >> (w-ptgt);
          len -= w;
        }
     if (len + ptgt > w) *tgt = old;
   }

