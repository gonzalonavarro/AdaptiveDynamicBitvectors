
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

	// supports leaf arrays of size up to 2^32-1

#include "leafId.h"

#define MaxBlockWords 128 // measured in uint64_t's
#define Gamma 0.75 // new blocks try to be this fraction full

       // size that a newly created leaf should have measured in elements

extern inline uint leafIdNewSize(uint width)

   { return leafIdMaxSize(width) * Gamma;
   }

        // max leaf size measured in elements

extern inline uint leafIdMaxSize(uint width)

   { return (MaxBlockWords*w/width / 2) * 2; // make it even
   }

	// creates an empty leafId of elements of width width

leafId leafIdCreate (uint width)

   { leafId B = (leafId)myalloc(sizeof(struct s_leafId));
     B->size = 0;
     B->width = width;
     B->isStat = 0;
     B->data = (uint64_t*)myalloc(MaxBlockWords * sizeof(uint64_t));
     return B;
   }

        // creates a static leafId of n elements of width width
        // from data in an already packed array, which is simply pointed

leafId leafIdCreateStaticFromPacked (uint64_t *data, uint n, uint width)

   { leafId B;
     B = (leafId)myalloc(sizeof(struct s_leafId));
     B->size = n;
     B->width = width;
     B->isStat = 1;
     B->data = data;
     return B;
   }

        // creates a dynamic leafId of n elements of width width
        // from data[i..] in an already packed array, which is not destroyed

leafId leafIdCreateFromPacked (uint64_t *data, uint64_t i, uint n, uint width)

   { leafId B;
     if (n == 0) return leafIdCreate(width);
     B = (leafId)myalloc(sizeof(struct s_leafId));
     B->size = n;
     B->width = width;
     B->isStat = 0;
     B->data = (uint64_t*)myalloc(MaxBlockWords * sizeof(uint64_t));
     if (i == 0) memcpy(B->data,data,(n*width+7)/8);
     else copyBits(B->data,0,data,i*width,n*width);
     return B;
   }

	// converts an array of n uint64_t into a leafId of n elements
	// of width width. data is pointed to and will be freed 

leafId leafIdCreateFrom64 (uint64_t *data, uint n, uint width, uint isStat)

   { uint i,p;
     leafId B;
     uint64_t word;
     if (n == 0) return leafIdCreate(width); // does not accept empty static
     B = (leafId)myalloc(sizeof(struct s_leafId));
     B->size = n;
     B->width = width;
     B->isStat = isStat;
     if (isStat)
	{ if (B->width == w) { B->data = data; return B; }
	  B->data = (uint64_t*)mycalloc(((n*width+w-1)/w),sizeof(uint64_t));
	}
     else B->data = (uint64_t*)mycalloc(MaxBlockWords,sizeof(uint64_t));
     p = 0;
     for (i=0;i<n;i++)
	 { word = data[i];
	   B->data[p/w] |= word << (p%w);
	   if ((p%w)+width > w) B->data[p/w+1] = word >> (w-(p%w));
	   p += width;
	 }
     free(data);
     return B;
   }

	// converts an array of n uint32_t into a leafId of n elements
	// of width width. data is pointed to and will be freed 

leafId leafIdCreateFrom32 (uint32_t *data, uint n, uint width, uint isStat)

   { uint i,p;
     leafId B;
     uint32_t word;
     if (n == 0) return leafIdCreate(width); // does not allow empty static
     B = (leafId)myalloc(sizeof(struct s_leafId));
     B->size = n;
     B->width = width;
     B->isStat = isStat;
     if (isStat)
	  B->data = (uint64_t*)mycalloc(((n*width+w-1)/w),sizeof(uint64_t));
     else B->data = (uint64_t*)mycalloc(MaxBlockWords,sizeof(uint64_t));
     p = 0;
     for (i=0;i<n;i++)
	 { word = data[i];
	   B->data[p/w] |= word << (p%w);
	   if ((p%w)+width > w) B->data[p/w+1] = word >> (w-(p%w));
	   p += width;
	 }
     free(data);
     return B;
   }

	// destroys B, frees data 

void leafIdDestroy (leafId B)

   { myfree(B->data);
     myfree(B);
   }

	// saves leaf data to file, which must be opened for writing

void leafIdSave (leafId B, FILE *file)

   { myfwrite(&B->size,sizeof(uint64_t),1,file);
     myfwrite(&B->width,sizeof(byte),1,file);
     myfwrite(&B->isStat,sizeof(byte),1,file);
     if (B->size != 0)
        myfwrite (B->data,sizeof(uint64_t),(B->size*B->width+w-1)/w,file);
   }

        // loads leaf data from file, which must be opened for reading
        // size is the number of elements

leafId leafIdLoad (FILE *file)

    { leafId B;
      uint64_t *data;
      uint64_t size;
      byte width,isStat;
      myfread(&size,sizeof(uint64_t),1,file);
      myfread(&width,sizeof(byte),1,file);
      myfread(&isStat,sizeof(byte),1,file);
      data = (uint64_t*)myalloc(((size*width+w-1)/w)*sizeof(uint64_t));
      myfread (data,sizeof(uint64_t),(size*width+w-1)/w,file);
      if (isStat) B = leafIdCreateStaticFromPacked(data,size,width);
      else { B = leafIdCreateFromPacked(data,0,size,width);
             myfree(data);
	   }
      return B;
    }

	// gives (allocated) space of B in w-bit words

uint leafIdSpace (leafId B)

   { return (sizeof(struct s_leafId)+sizeof(uint64_t)-1)/sizeof(uint64_t) 
	    + (B->isStat ? ((B->size*B->width+w-1)/w) : MaxBlockWords);
   }

	// gives array length

extern inline uint leafIdLength (leafId B)

   { return B->size;
   }

        // sets value for B[i]=v
	// assumes not static, i is right, and value fits in width

void leafIdWrite (leafId B, uint i, uint64_t v)

   { uint64_t word;
     uint iq,ir;
     if (B->width == w) { B->data[i] = v; return; }
     word = (((uint64_t)1) << B->width) - 1;
     iq = i*B->width/w;
     ir = (i*B->width)%w;
     B->data[iq] &= ~(word << ir);
     B->data[iq] |= v << ir;
     if (ir+B->width > w)
	{ B->data[iq+1] &= ~(word >> (w-ir));
	  B->data[iq+1] |= v >> (w-ir);
	}
   }

        // inserts v at B[i], assumes that insertion is possible
	// assumes not static, i is right, and value fits in width

void leafIdInsert (leafId B, uint i, uint64_t v)

   { uint nb = ++B->size*B->width/w;
     uint ib = i*B->width/w;
     uint ir = (i*B->width)%w;
     int b;

     if (B->width == w)
	{ for (b=nb;b>ib;b--) B->data[b] = B->data[b-1];
	  B->data[ib] = v;
	  return;
	}
     for (b=nb;b>ib;b--)
	 B->data[b] = (B->data[b] << B->width) | 
		      (B->data[b-1] >> (w-B->width));
     if (ir+B->width < w)
        B->data[ib] = (B->data[ib] & ((((uint64_t)1) << ir) - 1)) | (v << ir) |
                      ((B->data[ib] << B->width) & 
		       ((~(uint64_t)0) << (ir+B->width)));
     else {
	B->data[ib] = (B->data[ib] & ((((uint64_t)1) << ir) - 1)) | (v << ir);
	B->data[ib+1] = (B->data[ib+1] & ((~(uint64_t)0) << (ir+B->width-w))) |
			(v >> (w-ir));
        }
   }

        // deletes B[i]
	// assumes not static, i is right, and value fits in width

void leafIdDelete (leafId B, uint i)

   { uint nb = B->size--*B->width/w;
     uint ib = i*B->width/w;
     uint ir = (i*B->width)%w;
     int b;

     if (B->width == w)
	{ for (b=ib+1;b<nb;b++) B->data[b-1] = B->data[b];
	  return;
	}
     if (ir+B->width <= w)
        B->data[ib] = (B->data[ib] & ((((uint64_t)1) << ir) - 1)) |
		      ((B->data[ib] >> B->width) & (~((uint64_t)0) << ir));
     else {
	B->data[ib+1] = (B->data[ib+1] & (~(uint64_t)0) << (ir+B->width-w)) |
			((B->data[ib] >> (w-B->width)) & 
			 ((((uint64_t)1) << (ir+B->width-1)) - 1));
	B->data[ib] &= (~(uint64_t)0) << (w-B->width);
        }
     for (b=ib+1;b<=nb;b++) {
	 B->data[b-1] |= B->data[b] << (w-B->width);
	 B->data[b] >>= B->width;
	 }
   }

	// access B[i], assumes i is right

extern inline uint64_t leafIdAccess (leafId B, uint i)

   { uint64_t word;
     if (B->width == w) return B->data[i];
     uint iq = i*B->width/w;
     uint ir = (i*B->width)%w;
     word = (B->data[iq] >> ir);
     if (ir+B->width > w) word |= B->data[iq+1] << (w-ir);
     return word & ((((uint64_t)1) << B->width) - 1);
   }

        // read bits [i..i+l-1], onto D[0...] of 64 bits

void leafIdRead64 (leafId B, uint i, uint l, uint64_t *D)

   { uint j,width,iq,ir;
     width = B->width;
     if (width == w) 
	{ memcpy(D,B->data+i,l*sizeof(uint64_t));
	  return;
	}
     iq = i*B->width/w;
     ir = (i*B->width)%w;
     for (j=0;j<l;j++)
	{ D[j] = B->data[iq] >> ir;
	  if (ir+width >= w) 
	     { iq++;
	       D[j] |= B->data[iq] << (w-ir);
	       ir += width-w;
	     }
	  else ir += width; 
	  D[j] &= ((((uint64_t)1) << width) - 1);
	}
   }

        // read bits [i..i+l-1], onto D[0...] of 32 bits

void leafIdRead32 (leafId B, uint i, uint l, uint32_t *D)

   { uint j,width,iq,ir;
     width = B->width;
     iq = i*B->width/w;
     ir = (i*B->width)%w;
     for (j=0;j<l;j++)
	{ D[j] = B->data[iq] >> ir;
	  if (ir+width >= w) 
	     { iq++;
	       D[j] |= B->data[iq] << (w-ir);
	       ir += width-w;
	     }
	  else ir += width; 
	  D[j] &= ((((uint64_t)1) << width) - 1);
	}
   }

        // finds first value >= c in [i..j], which must be increasing
        // returns j+1 if not found

uint leafIdNext (leafId B, uint i, uint j, uint64_t c)

   { uint width,iq,ir,d,m;
     uint64_t v;
     width = B->width;
	// answer is i
     iq = i*B->width/w;
     ir = (i*B->width)%w;
     v = B->data[iq] >> ir;
     if (ir+width >= w) 
	v |= B->data[iq+1] << (w-ir);
     v &= ((((uint64_t)1) << width) - 1);
     if (v >= c) return i;
	// answer is j+1
     iq = j*B->width/w;
     ir = (j*B->width)%w;
     v = B->data[iq] >> ir;
     if (ir+width >= w) 
	v |= B->data[iq+1] << (w-ir);
     v &= ((((uint64_t)1) << width) - 1);
     if (v < c) return j+1;
	// invariant data[i] < c and data[j] >= c
     d = 1;
     while (i+d <= j)
        { iq = (i+d)*B->width/w;
          ir = ((i+d)*B->width)%w;
	  v = B->data[iq] >> ir;
	  if (ir+width >= w) 
	     v |= B->data[iq+1] << (w-ir);
	  v &= ((((uint64_t)1) << width) - 1);
	  if (v >= c) break;
          i += d; d <<= 1; 
	}
     d = min(j,i+d); // data[d] >= j
     while (i+1<d)
        { m = (i+d)>>1;
          iq = m*B->width/w;
          ir = (m*B->width)%w;
	  v = B->data[iq] >> ir;
	  if (ir+width >= w) 
	     v |= B->data[iq+1] << (w-ir);
	  v &= ((((uint64_t)1) << width) - 1);
          if (v < c) i = m; else d = m;
        }
     return d;
   }


