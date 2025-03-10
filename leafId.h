
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

#ifndef INCLUDEDleafId
#define INCLUDEDleafId

	// supports leaf arrays of size up to 2^32-1

#include "basics.h"

typedef struct s_leafId
   { uint64_t size; // elements represented
     byte width; // bits used per element, up to w
     byte isStat; // does not accept indels
     uint64_t *data; 
   } *leafId;

	// size that a newly created leaf should have, and max leaf size
	// measured in elements
extern inline uint leafIdNewSize(uint width);
extern inline uint leafIdMaxSize(uint width);

	// creates an empty leafId with elements of width width
leafId leafIdCreate (uint width);

        // creates a static leafId of n elements of width width
        // from data in an already packed array, which is simply pointed
leafId leafIdCreateStaticFromPacked (uint64_t *data, uint n, uint width);

        // creates a dynamic leafId of n elements of width width
        // from data[i..] in an already packed array, which is not destroyed
leafId leafIdCreateFromPacked (uint64_t *data, uint64_t i, uint n, uint width);

	// converts an array of n uint64_t into a leafId of n elements
	// of width width. data is pointed to and will be freed 
leafId leafIdCreateFrom64 (uint64_t *data, uint n, uint width, uint isStat);

	// converts an array of n uint32_t into a leafId of n elements
	// of width width. data is pointed to and will be freed 
leafId leafIdCreateFrom32 (uint32_t *data, uint n, uint width, uint isStat);

	// destroys B, frees data 
void leafIdDestroy (leafId B);

	// saves leaf data to file, which must be opened for writing
void leafIdSave (leafId B, FILE *file);

        // loads leaf data from file, which must be opened for reading
leafId leafIdLoad (FILE *file);

	// gives (allocated) space of B in w-bit words
uint leafIdSpace (leafId B);

	// gives array length
extern inline uint leafIdLength (leafId B);

        // sets value for B[i] = v
	// assumes not static, i is right, and value fits in width
void leafIdWrite (leafId B, uint i, uint64_t v);

        // inserts v at B[i], assumes that insertion is possible
	// assumes not static, i is right, and value fits in width
void leafIdInsert (leafId B, uint i, uint64_t v);

        // deletes B[i], assumes i is right
	// assumes not static, i is right, and value fits in width
void leafIdDelete (leafId B, uint i);

	// access B[i], assumes i is right
extern inline uint64_t leafIdAccess (leafId B, uint i);

        // read bits [i..i+l-1], onto D[0...] of 64 bits
void leafIdRead64 (leafId B, uint i, uint l, uint64_t *D);

        // read bits [i..i+l-1], onto D[0...] of 32 bits
void leafIdRead32 (leafId B, uint i, uint l, uint32_t *D);

	// finds first value >= c in [i..j], which must be increasing
	// returns j+1 if not found
uint leafIdNext (leafId B, uint i, uint j, uint64_t c);

#endif
