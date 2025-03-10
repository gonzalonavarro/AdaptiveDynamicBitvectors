
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

#ifndef INCLUDEDleafBV
#define INCLUDEDleafBV

	// supports leaf bitvectors of size up to 2^32-1

#include "basics.h"

typedef struct s_leafBV
   { uint size; // bits represented
     uint ones; // # 1s
     uint64_t *data; 
   } *leafBV;

	// size that a newly created leaf should have, and max leaf size
	// multiples of w
extern inline uint leafNewSize(void);
extern inline uint leafMaxSize(void);

	// creates an empty leafBV 
leafBV leafCreate (void);

	// converts a bit array into a leafBV of n bits
	// data is copied but freed only if freeit
leafBV leafCreateFrom (uint64_t *data, uint n, int freeit);

	// destroys B, frees data 
void leafDestroy (leafBV B);

	// saves leaf data to file, which must be opened for writing
void leafSave (leafBV B, FILE *file);

        // loads leaf data from file, which must be opened for reading
        // size is the number of bits
leafBV leafLoad (FILE *file, uint size);

	// gives (allocated) space of B in w-bit words
uint leafSpace (leafBV B);

	// gives bit length
extern inline uint leafLength (leafBV B);

	// gives number of 1s
extern inline uint leafOnes (leafBV B);

       // sets value for B[i]= (v != 0), assumes i is right
        // returns difference in 1s
int leafWrite (leafBV B, uint i, uint v);

        // inserts v at B[i], assumes i is right and that insertion is possible
void leafInsert (leafBV B, uint i, uint v);

        // deletes B[i], assumes i is right
        // returns difference in 1s
int leafDelete (leafBV B, uint i);

	// access B[i], assumes i is right
extern inline uint leafAccess (leafBV B, uint i);

        // read bits [i..i+l-1], onto D[j...]
void leafRead (leafBV B, uint i, uint l, uint64_t *D, uint64_t j);

	// computes rank_1(B,i), zero-based, assumes i is right
uint leafRank (leafBV B, uint i);

	// computes select_1(B,j), zero-based, assumes j is right
uint leafSelect (leafBV B, uint j);

	// computes select_0(B,j), zero-based, assumes j is right
uint leafSelect0 (leafBV B, uint j);

        // computes next_1(B,i), zero-based and including i
        // returns -1 if no answer

int leafNext (leafBV B, uint i);

        // computes next_0(B,i), zero-based and including i
        // returns -1 if no answer

int leafNext0 (leafBV B, uint i);

#endif
