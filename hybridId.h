
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

#ifndef INCLUDEDhybridId
#define INCLUDEDhybridId

	// supports hybrid arrays of size up to 2^64-1

#include "leafId.h"
#include "hybridBV.h" // to include nodeType

typedef struct s_hybridId *hybridId;

typedef struct s_dynamicId
   { uint64_t size;
     byte width; // up to w
     uint64_t leaves; // leaves below node
     uint64_t accesses; // since last update
     hybridId left,right; // hybridIds
   } *dynamicId;

typedef struct s_hybridId
   { nodeType type;
     union
      { leafId stat;
        leafId leaf;
        dynamicId dyn;
      } bv;
   } *hybridId;
      
extern float ThetaId; // reconstruction factor

	// creates an empty hybridId, of width width
hybridId hybridIdCreate (uint width);

	// converts an array of uint64_t into a hybridId of n elements
	// of width width. data is pointed to and will be freed 
hybridId hybridIdCreateFrom64 (uint64_t *data, uint64_t n, uint width);

	// converts an array of uint32_t into a hybridId of n elements
	// of width width. data is pointed to and will be freed 
hybridId hybridIdCreateFrom32 (uint32_t *data, uint64_t n, uint width);

	// destroys B, frees data 
void hybridIdDestroy (hybridId B);

	// writes B to file, which must be opened for writing
void hybridIdSave (hybridId B, FILE *file);

	// loads hybridId from file, which must be opened for reading
hybridId hybridIdLoad (FILE *file);

	// gives space of hybridId in w-bit words
uint64_t hybridIdSpace (hybridId B);

	// gives number of elements length
extern inline uint64_t hybridIdLength (hybridId B);

	// gives width of elements 
extern inline uint hybridIdWidth (hybridId B);

	// sets value for B[i] = v, assumes i is right and v fits in width
void hybridIdWrite (hybridId B, uint64_t i, uint64_t v);

	// inserts v at B[i], assumes i is right and v fits in width
void hybridIdInsert (hybridId B, uint64_t i, uint64_t v);

	// deletes B[i], assumes i is right
void hybridIdDelete (hybridId B, uint64_t i);

	// access B[i], assumes i is right
uint64_t hybridIdAccess (hybridId B, uint64_t i);

        // read values [i..i+l-1], onto D[0...], of uint64_t
void hybridIdRead64 (hybridId B, uint64_t i, uint64_t l, uint64_t *D);

        // read values [i..i+l-1], onto D[0...], of uint32_t
void hybridIdRead32 (hybridId B, uint64_t i, uint64_t l, uint32_t *D);

#endif
