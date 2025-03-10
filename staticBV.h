
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

#ifndef INCLUDEDstaticBV
#define INCLUDEDstaticBV

	// supports static bitvectors of size up to 2^64-1

#include "basics.h"

typedef struct s_staticBV {
    uint64_t size; // number of bits
    uint64_t ones; // number of 1s
    uint64_t* data; // the bits
    uint64_t *S; // superblocks
    uint16_t *B; // blocks
    } *staticBV;

	// converts a bit array into a bitvector of n bits
	// data is pointed to and will be freed 
staticBV staticCreateFrom (uint64_t *data, uint64_t n);

	// destroys B, frees data 
void staticDestroy (staticBV B);

	// writes B's data to file, which must be opened for writing
void staticSave (staticBV B, FILE *file);

	// loads bitvector's data from file, which must be opened for reading
	// size is the number of bits
staticBV staticLoad (FILE *file, uint64_t size);

	// gives space of bitvector in w-bit words
uint64_t staticSpace (staticBV B);

        // data of staticBV
extern inline uint64_t *staticBits (staticBV B);

	// gives bit length
extern inline uint64_t staticLength (staticBV B);

	// gives number of ones
extern inline uint64_t staticOnes (staticBV B);

	// access B[i], assumes i is right
extern inline uint staticAccess (staticBV B, uint64_t i);

        // read bits [i..i+l-1], writes onto D[j..]
extern inline void staticRead (staticBV B, uint64_t i, uint64_t l, 
			       uint64_t *D, uint64_t j);

	// computes rank_1(B,i), zero-based, assumes i is right
extern inline uint64_t staticRank (staticBV B, uint64_t i);

	// computes rank_0(B,i), zero-based, assumes i is right
extern inline uint64_t staticRank0 (staticBV B, uint64_t i);

	// computes select_1(B,j), zero-based, assumes j is right
extern uint64_t staticSelect (staticBV B, uint64_t j);

	// computes select_0(B,j), zero-based, assumes j is right
extern uint64_t staticSelect0 (staticBV B, uint64_t j);

	// computes next_1(B,i), zero-based, assumes i is right
extern int64_t staticNext (staticBV B, uint64_t i);

	// computes next_0(B,i), zero-based, assumes i is right
extern int64_t staticNext0 (staticBV B, uint64_t i);

#endif
