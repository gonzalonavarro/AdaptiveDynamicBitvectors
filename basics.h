
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

#ifndef INCLUDEDbasics
#define INCLUDEDbasics

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned char byte;

typedef uint32_t uint;

#define w (8*sizeof(uint64_t))
#define min(x,y) (((x)<(y))?(x):(y))
#define max(x,y) (((x)>(y))?(x):(y))

void *myalloc (size_t n);
void *mycalloc (size_t n, size_t s);
void *myrealloc (void *p, size_t n);
void myfree (void *p);

void myfread(void *ptr, size_t size, size_t nmemb, FILE *stream);
void myfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);

	// number of bits needed to represent n, gives 1 for n=0
uint numbits (uint n);

        // counts # of 1s in y
extern inline uint popcount (uint64_t y);

        // copies len bits starting at *src + psrc
        // to tgt from bit position ptgt
        // WARNING: leave at least one extra word to spare in tgt
void copyBits (uint64_t *tgt, uint64_t ptgt,
               uint64_t *src, uint64_t psrc, uint64_t len);

#endif
