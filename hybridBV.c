
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

	// supports hybrid bitvectors of size up to 2^64-1

#include "hybridBV.h"

extern float Theta = 0.01; // Theta * length reads => rebuild as static

static const float Epsilon = 0.1; // do not flatten leaves of size over Epsilon * n

static const float Alpha = 0.65; // balance factor 3/5 < . < 1

static const float TrfFactor = 0.125; // TrfFactor * MaxLeafSize to justify transferLeft/Right

static const int MinLeavesToBalance = 5; // min number of leaves to balance the tree

static const float MinFillFactor = 0.3; // less than this involves rebuild. 
				// Must be <= Gamma/2

	// internal, to study behavior
extern uint64_t flattenMax = 0;
extern uint64_t flattenAccess = 0;
extern uint64_t flattenBalance = 0;
extern uint64_t flattenFill = 0;

static inline mustFlatten(hybridBV B, uint64_t n)
  { return ((B->bv.dyn->size <= Epsilon * n) && (B->bv.dyn->accesses >= Theta * B->bv.dyn->size));
    // return (B->bv.dyn->accesses >= Theta * B->bv.dyn->size);
  }

	// creates an empty hybridBV

hybridBV hybridCreate (void)

   { hybridBV B = myalloc(sizeof(struct s_hybridBV));
     B->type = tLeaf;
     B->bv.leaf = leafCreate();
     return B;
   }

	// converts a bit array into a hybridBV of n bits
	// data is pointed to and will be freed 

hybridBV hybridCreateFrom (uint64_t *data, uint64_t n)

   { hybridBV B = myalloc(sizeof(struct s_hybridBV));
     if (n > leafNewSize()*w)
	{ B->type = tStatic;
          B->bv.stat = staticCreateFrom(data,n);
	}
     else
	{ B->type = tLeaf;
          B->bv.leaf = leafCreateFrom(data,n,1);
	}
     return B;
   }

	// destroys B, frees data 

void hybridDestroy (hybridBV B)

   { if (B->type == tLeaf) leafDestroy(B->bv.leaf);
     else if (B->type == tStatic) staticDestroy(B->bv.stat);
     else { hybridDestroy(B->bv.dyn->left);
            hybridDestroy(B->bv.dyn->right);
	    myfree(B->bv.dyn);
	  }
     myfree(B);
   }

	// creates a static version of B, rewriting it but not its address

	// version of hybridRead that does not count accesses, for internal use

static void myread (hybridBV B, uint64_t i, uint64_t l, uint64_t *D, uint64_t j)

   { uint64_t lsize;
     if (B->type == tLeaf) { leafRead(B->bv.leaf,i,l,D,j); return; }
     if (B->type == tDynamic)
        { lsize = hybridLength(B->bv.dyn->left);
     	  if (i+l < lsize) myread(B->bv.dyn->left,i,l,D,j);
	  else if (i >= lsize) myread(B->bv.dyn->right,i-lsize,l,D,j);
	  else { myread(B->bv.dyn->left,i,lsize-i,D,j);
		 myread(B->bv.dyn->right,0,l-(lsize-i),D,j+(lsize-i));
	       }
	  return;
          }
     staticRead(B->bv.stat,i,l,D,j);
   }

	// collects all the descending bits into an array, destroys bv.dyn

static uint64_t* collect (hybridBV B, uint64_t len)

   { uint64_t *D;
     D = (uint64_t*)myalloc(((len+w-1)/w)*sizeof(uint64_t));
     myread (B,0,len,D,0);
     hybridDestroy(B->bv.dyn->left);
     hybridDestroy(B->bv.dyn->right);
     myfree(B->bv.dyn);
     return D;
   }

	// gives number of leaves

inline uint64_t hybridLeaves (hybridBV B)

   { if (B->type == tLeaf) return 1;
     if (B->type == tStatic) return (staticLength(B->bv.stat)+leafNewSize()*w-1)
					/ (leafNewSize()*w);
     return B->bv.dyn->leaves;
   }

	// converts into a leaf if it's short or into a static otherwise
	// delta gives the difference in leaves (new - old)

static void flatten (hybridBV B, int64_t *delta)

   { uint64_t len;
     uint64_t *D;
     if (B->type != tDynamic) return;
     len = hybridLength(B);
     flattenAccess += len;
     if (len > flattenMax) flattenMax = len;
     *delta = - hybridLeaves(B);
     D = collect(B,len);
     if (len > leafNewSize()*w) // creates a static
        { B->type = tStatic;
          B->bv.stat = staticCreateFrom(D,len);
	}
     else
        { B->type = tLeaf;
          B->bv.leaf = leafCreateFrom(D,len,1);
	}
     *delta += hybridLeaves(B);
   }

	// splits a full leaf into two
	// returns a dynamicBV and destroys B

static dynamicBV splitLeaf (leafBV B)

   { hybridBV HB1,HB2;
     dynamicBV DB;
     leafBV LB;
     uint bsize; 

     bsize = (B->size/2+7)/8; // byte size of new left leaf
     LB = leafCreateFrom(B->data,bsize*8,0);
     HB1 = (hybridBV)myalloc(sizeof(struct s_hybridBV));
     HB1->type = tLeaf;
     HB1->bv.leaf = LB;
     LB = leafCreateFrom((uint64_t*)(((byte*)B->data)+bsize),B->size-bsize*8,0);
     HB2 = (hybridBV)myalloc(sizeof(struct s_hybridBV));
     HB2->type = tLeaf;
     HB2->bv.leaf = LB;
     DB = (dynamicBV)myalloc(sizeof(struct s_dynamicBV));
     DB->size = B->size;
     DB->ones = B->ones;
     DB->leaves = 2;
     DB->accesses = 0;
     DB->left = HB1;
     DB->right = HB2;
     leafDestroy(B);
     return DB;
   }

	// halves a static bitmap into leaves, leaving a leaf covering i
	// returns a dynamicBV and destroys B

static dynamicBV splitFrom (uint64_t *data, uint64_t n, uint64_t ones,
			    uint64_t i)

   { hybridBV HB;
     dynamicBV DB,finalDB;
     uint bsize; // byte size of blocks to create
     uint blen; // bit size of block to create
     uint64_t *segment;
     uint64_t nblock;
     byte *start,*mid,*end;

     HB = NULL;
     blen = leafNewSize() * w;
     bsize = (blen+7)/8; // +7 not really needed as 8 | w
     nblock = (n+blen-1)/blen; // total blocks 
     start = (byte*)data;
     end = start + (n+7)/8;
     while (nblock >= 2) {
        DB = (dynamicBV)myalloc(sizeof(struct s_dynamicBV));
	if (HB == NULL) finalDB = DB;
	else { 
           HB->type = tDynamic;
           HB->bv.dyn = DB;
	   }
        DB->size = n;
	DB->ones = ones;
        DB->leaves = nblock;
        DB->accesses = 0;
	mid = start+(nblock/2)*bsize;
     	if (i < (nblock/2)*blen) { // split the left half
		// create right half
           DB->right = HB = (hybridBV)myalloc(sizeof(struct s_hybridBV));
	   if (n - (nblock/2)*blen > leafNewSize() * w) { // create a static
	      segment = (uint64_t*)myalloc(end-mid);
	      memcpy(segment,mid,end-mid);
	      HB->type = tStatic;
	      HB->bv.stat = staticCreateFrom(segment,n-(nblock/2)*blen);
	      }
	   else { // create a leaf
	      HB->type = tLeaf;
	      HB->bv.leaf =leafCreateFrom((uint64_t*)mid,n-(nblock/2)*blen,0);
	      }
		// continue on left half
	   end = mid;
	   nblock = nblock/2;
	   n = nblock * blen;
	   ones -= hybridOnes(HB);
           DB->left = HB = (hybridBV)myalloc(sizeof(struct s_hybridBV));
	   }
	else { // split the right half
		// create left half
           DB->left = HB = (hybridBV)myalloc(sizeof(struct s_hybridBV));
	   if ((nblock/2)*blen > leafNewSize() * w) { // create a static
	      segment = (uint64_t*)myalloc(mid-start);
	      memcpy(segment,start,mid-start);
	      HB->type = tStatic;
	      HB->bv.stat = staticCreateFrom(segment,(nblock/2)*blen);
	      }
	   else { // create a leaf
	      HB->type = tLeaf;
	      HB->bv.leaf =leafCreateFrom((uint64_t*)start,(nblock/2)*blen,0);
	      }
		// continue for right half
	   start = mid;
	   n = n-(nblock/2)*blen;
	   i = i-(nblock/2)*blen;
	   ones -= hybridOnes(HB);
	   nblock = nblock - nblock/2;
           DB->right = HB = (hybridBV)myalloc(sizeof(struct s_hybridBV));
	   }
	}
	// finally, the leaf where i lies
     HB->type = tLeaf;
     HB->bv.leaf =leafCreateFrom((uint64_t*)start,n,0);
     return finalDB;
   }

static dynamicBV split (staticBV B, uint64_t i)

   { dynamicBV DB;
     DB = splitFrom(B->data,staticLength(B),staticOnes(B),i);
     staticDestroy(B);
     return DB;
   }

	// balance by rebuilding: flattening + splitting
	// assumes B is dynamic

static inline int canBalance (uint64_t n, int dleft, int dright)

   { uint b = leafNewSize() * w; // bit size of leaves to create in split
     uint64_t left = (((n+b-1)/b)/2)*b; // bit size of left part 
     uint64_t right = n-left; // bit size of right part 
//     if (left+dleft < (1-Alpha)*(n+dleft+dright)) return 0;
//     if (right+dright < (1-Alpha)*(n+dleft+dright)) return 0;
     if (left+dleft > Alpha*(n+dleft+dright)) return 0;
     if (right+dright > Alpha*(n+dleft+dright)) return 0;
     return 1;
   }

static void balance (hybridBV B, uint64_t i, int64_t *delta)

   { uint64_t len = hybridLength(B);
     uint64_t ones = hybridOnes(B);
     uint64_t *D;
     flattenBalance += len;
     *delta = - hybridLeaves(B);
     D = collect(B,len);
     B->bv.dyn = splitFrom(D,len,ones,i);
     *delta += hybridLeaves(B);
     myfree(D);
   }

	// merge the two leaf children of B into a leaf
	// returns a leafBV and destroys B

static leafBV mergeLeaves (dynamicBV B)

   { leafBV LB1,LB2;

     LB1 = B->left->bv.leaf;
     LB2 = B->right->bv.leaf;
     copyBits(LB1->data,LB1->size,LB2->data,0,LB2->size);
     LB1->size += LB2->size;
     LB1->ones += LB2->ones;
     leafDestroy(LB2);
     myfree(B);
     return LB1;
   }

	// transfers bits from the right to the left leaves children of B
	// tells if it transfered something

static int transferLeft (dynamicBV B)

   { leafBV LB1,LB2;
     uint i,trf,ones,words;
     uint64_t *segment;

     LB1 = B->left->bv.leaf;
     LB2 = B->right->bv.leaf;
     trf = (LB2->size-LB1->size+1)/2;
     if (trf < leafMaxSize() * w * TrfFactor) return 0;
     copyBits(LB1->data,LB1->size,LB2->data,0,trf);
     LB1->size += trf;
     LB2->size -= trf;
     segment = (uint64_t*)myalloc(leafMaxSize()*sizeof(uint64_t));
     copyBits(segment,0,LB2->data,trf,LB2->size);
     words = trf / w;
     ones = 0;
     for (i=0;i<words;i++) ones += popcount(LB2->data[i]);
     if (trf % w) ones += popcount(LB2->data[words] & 
				   (((uint64_t)1) << (trf % w)) - 1);
     LB1->ones += ones;
     LB2->ones -= ones;
     memcpy(LB2->data,segment,(LB2->size+7)/8);
     myfree(segment);
     return 1;
   }

	// transfers bits from the left to the right leaves children of B
	// tells if it transfered something

static int transferRight (dynamicBV B)

   { leafBV LB1,LB2;
     uint i,trf,ones,words;
     uint64_t *segment;

     LB1 = B->left->bv.leaf;
     LB2 = B->right->bv.leaf;
     trf = (LB1->size-LB2->size+1)/2;
     if (trf < leafMaxSize() * w * TrfFactor) return 0;
     segment = (uint64_t*)myalloc(leafMaxSize()*sizeof(uint64_t));
     memcpy(segment,LB2->data,(LB2->size+7)/8);
     copyBits(LB2->data,0,LB1->data,LB1->size-trf,trf);
     words = trf / w;
     ones = 0;
     for (i=0;i<words;i++) ones += popcount(LB2->data[i]);
     if (trf % w) ones += popcount(LB2->data[words] & 
				   (((uint64_t)1) << (trf % w)) - 1);
     LB1->ones -= ones;
     LB2->ones += ones;
     copyBits(LB2->data,trf,segment,0,LB2->size);
     LB1->size -= trf;
     LB2->size += trf;
     myfree(segment);
     return 1;
   }

	// writes B to file, which must be opened for writing

void hybridSave (hybridBV B, FILE *file)

   { int64_t delta;
     flatten(B,&delta);
     if (B->type == tStatic) 
	{ myfwrite (&B->bv.stat->size,sizeof(uint64_t),1,file);
	  staticSave(B->bv.stat,file);
	}
     else 
	{ myfwrite (&B->bv.leaf->size,sizeof(uint64_t),1,file);
	  leafSave(B->bv.leaf,file);
	}
   }

	// loads hybridBV from file, which must be opened for reading

hybridBV hybridLoad (FILE *file)

   { uint64_t size;
     hybridBV B = myalloc(sizeof(struct s_hybridBV));
     myfread (&size,sizeof(uint64_t),1,file);
     if (size > leafNewSize()*w)
        { B->type = tStatic;
          B->bv.stat = staticLoad(file,size);
	}
     else
        { B->type = tLeaf;
          B->bv.leaf = leafLoad(file,size);
	}
     return B;
   }

	// gives space of hybridBV in w-bit words

uint64_t hybridSpace (hybridBV B)

   { uint64_t s;
     s = (sizeof(struct s_hybridBV)*8+w-1)/w;
     if (B->type == tLeaf) return s+leafSpace(B->bv.leaf);
     else if (B->type == tStatic) return s+staticSpace(B->bv.stat);
     return s+(sizeof(struct s_dynamicBV)*8+w-1)/w+
	    hybridSpace(B->bv.dyn->left)+hybridSpace(B->bv.dyn->right);
   }

	// gives bit length

inline uint64_t hybridLength (hybridBV B)

   { if (B->type == tLeaf) return leafLength(B->bv.leaf);
     if (B->type == tStatic) return staticLength(B->bv.stat);
     return B->bv.dyn->size;
   }

	// gives number of ones

inline uint64_t hybridOnes (hybridBV B)

   { if (B->type == tLeaf) return leafOnes(B->bv.leaf);
     if (B->type == tStatic) return staticOnes(B->bv.stat);
     return B->bv.dyn->ones;
   }

	// sets value for B[i]= (v != 0), assumes i is right
	// returns the difference in 1s

int hybridWrite (hybridBV B, uint64_t i, uint v)

   { uint64_t lsize;
     int dif;
     if (B->type == tStatic) { 
	B->type = tDynamic;
	B->bv.dyn = split(B->bv.stat,i); // does not change #leaves!
	}
     if (B->type == tLeaf) 
	return leafWrite(B->bv.leaf,i,v);
     B->bv.dyn->accesses = 0; // reset
     lsize = hybridLength(B->bv.dyn->left);
     if (i < lsize) dif = hybridWrite(B->bv.dyn->left,i,v);
     else dif = hybridWrite(B->bv.dyn->right,i-lsize,v);
     B->bv.dyn->ones += dif;
     return dif;
   }

	// changing leaves is uncommon and only then we need to recompute
	// leaves. we do our best to avoid this overhead in typical operations

static void irecompute (hybridBV B, uint64_t i)

   { uint64_t lsize;
     if (B->type == tDynamic) {
        lsize = hybridLength(B->bv.dyn->left);
        if (i < lsize) irecompute(B->bv.dyn->left,i);
        else irecompute(B->bv.dyn->right,i-lsize);
        B->bv.dyn->leaves = hybridLeaves(B->bv.dyn->left) +
			    hybridLeaves(B->bv.dyn->right);
	}
   }

static void rrecompute (hybridBV B, uint64_t i, uint64_t l)

   { uint64_t lsize;
     if (B->type == tDynamic) {
	lsize = hybridLength(B->bv.dyn->left);
     	if (i+l < lsize) rrecompute(B->bv.dyn->left,i,l);
	else if (i >= lsize) rrecompute(B->bv.dyn->right,i-lsize,l);
	else { rrecompute(B->bv.dyn->left,i,lsize-i);
	       rrecompute(B->bv.dyn->right,0,l-(lsize-i));
	     }
        B->bv.dyn->leaves = hybridLeaves(B->bv.dyn->left) +
			    hybridLeaves(B->bv.dyn->right);
	}
   }

	// inserts v at B[i], assumes i is right

static void insert (hybridBV B, uint64_t i, uint v, uint *recalc)

   { uint64_t lsize,rsize;
     int64_t delta;
     if (B->type == tStatic) { 
	B->type = tDynamic;
	B->bv.dyn = split(B->bv.stat,i); // does not change #leaves!
	}
     if (B->type == tLeaf) {
	if (leafLength(B->bv.leaf) == leafMaxSize() * w) // split
	   { B->type = tDynamic;
	     B->bv.dyn = splitLeaf(B->bv.leaf);
	     *recalc = 1; // leaf added
	   }
	else 
	   { leafInsert(B->bv.leaf,i,v);
	     return;
	   }
	}
     B->bv.dyn->accesses = 0; // reset
     lsize = hybridLength(B->bv.dyn->left);
     rsize = hybridLength(B->bv.dyn->right);
     if (i < lsize) {  // insert on left child
        if ((lsize == leafMaxSize() * w)     // will overflow if leaf
	     && (rsize < leafMaxSize() * w)  // can avoid if leaf
	     && (B->bv.dyn->left->type == tLeaf) // both are leaves
	     && (B->bv.dyn->right->type == tLeaf) 
	     && transferRight(B->bv.dyn)) { // avoided, transferred to right
	   insert(B,i,v,recalc); // now could be to the right!
	   return;
	   }
	if ((lsize+1 > Alpha*(lsize+rsize+1))
	    && (lsize+rsize >= MinLeavesToBalance*leafMaxSize()*w) 
	    && canBalance(lsize+rsize,1,0)) { // too biased
	   delta = 0;
	   balance(B,i,&delta);
	   if (delta) *recalc = 1;
	   insert(B,i,v,recalc); 
	   return;
	   }
	insert(B->bv.dyn->left,i,v,recalc); // normal recursive call 
	}
     else { // insert on right child
	if ((rsize == leafMaxSize() * w)    // will overflow if leaf
	     && (lsize < leafMaxSize() * w) // can avoid if leaf
	     && (B->bv.dyn->left->type == tLeaf) // both are leaves
	     && (B->bv.dyn->right->type == tLeaf) 
	     && transferLeft(B->bv.dyn)) { // avoided, transferred to right
	   insert(B,i,v,recalc); // now could be to the left!
	   return;
	   }
	if ((rsize+1 > Alpha*(lsize+rsize+1))
	    && (lsize+rsize >= MinLeavesToBalance*leafMaxSize()*w) 
	    && canBalance(lsize+rsize,0,1)) { // too biased
	   delta = 0;
	   balance(B,i,&delta);
	   if (delta) *recalc = 1;
	   insert(B,i,v,recalc);
	   return;
	   }
        insert(B->bv.dyn->right,i-lsize,v,recalc); // normal rec call
	}
     B->bv.dyn->size++;
     B->bv.dyn->ones += v;
   }

void hybridInsert (hybridBV B, uint64_t i, uint v)

   { uint recalc = 0;
     insert(B,i,v,&recalc);
     if (recalc) irecompute(B,i); // we went to the leaf now holding i
   }

	// deletes B[i], assumes i is right
	// returns difference in 1s

static int delete (hybridBV B, uint64_t i, uint *recalc)

   { uint64_t lsize,rsize;
     hybridBV B2;
     int dif;
     int64_t delta;
     if (B->type == tStatic) { 
	B->type = tDynamic;
	B->bv.dyn = split(B->bv.stat,i); // does not change #leaves!
	}
     if (B->type == tLeaf) 
	return leafDelete(B->bv.leaf,i);
     B->bv.dyn->accesses = 0; // reset
     lsize = hybridLength(B->bv.dyn->left);
     rsize = hybridLength(B->bv.dyn->right);
     if (i < lsize) { 
        if ((rsize > Alpha*(lsize+rsize-1))
	    && (lsize+rsize >= MinLeavesToBalance*leafMaxSize()*w) 
	    && canBalance(lsize+rsize,-1,0)) { // too biased
	   delta = 0;
	   balance(B,i,&delta); 
	   if (delta) *recalc = 1;
	   return delete(B,i,recalc); // now could enter in the right child!
	   }
	dif = delete(B->bv.dyn->left,i,recalc); // normal recursive call otherw
        if (lsize == 1) { // left child is now of size zero, remove
           hybridDestroy(B->bv.dyn->left);
           B2 = B->bv.dyn->right;
           myfree(B->bv.dyn);
           *B = *B2;
	   *recalc = 1; 
           return dif;
           } 
	}
     else {
        if ((lsize > Alpha*(lsize+rsize-1))
	    && (lsize+rsize >= MinLeavesToBalance*leafMaxSize()*w) 
	    && canBalance(lsize+rsize,0,-1)) { // too biased
	   delta = 0;
	   balance(B,i,&delta);
	   if (delta) *recalc = 1; 
	   return delete(B,i,recalc); // now could enter in the left child!
	   }
        dif = delete(B->bv.dyn->right,i-lsize,recalc); // normal recursive call
        if (rsize == 1) { // right child now size zero, remove
           hybridDestroy(B->bv.dyn->right);
           B2 = B->bv.dyn->left;
           myfree(B->bv.dyn);
           *B = *B2;
	   *recalc = 1; 
           return dif;
           }
	}
     B->bv.dyn->size--;
     B->bv.dyn->ones += dif;
     if (B->bv.dyn->size <= leafNewSize() * w) { // merge, must be leaves
	B->bv.leaf = mergeLeaves(B->bv.dyn);
	B->type = tLeaf;
	*recalc = 1; 
	}
     else if (B->bv.dyn->size < 
	      B->bv.dyn->leaves * leafNewSize() * w * MinFillFactor) {
	delta = 0;
	flattenFill += B->bv.dyn->size;
	flattenAccess -= B->bv.dyn->size;
	flatten(B,&delta); 
	if (delta) *recalc = 1;
	}
     return dif;
   }

int hybridDelete (hybridBV B, uint64_t i)

   { uint recalc = 0;
     int dif = delete(B,i,&recalc);
     if (recalc) { // the node is now at i-1 or at i, hard to know
        irecompute(B,i-1);
        irecompute(B,i); 
	}
     return dif;
   }

	// flattening is uncommon and only then we need to recompute
	// leaves. we do our best to avoid this overhead in typical queries

static void recompute (hybridBV B, uint64_t i, int64_t delta)

   { uint64_t lsize;
     if (B->type == tDynamic) {
        B->bv.dyn->leaves += delta;
        lsize = hybridLength(B->bv.dyn->left);
        if (i < lsize) recompute(B->bv.dyn->left,i,delta);
        else recompute(B->bv.dyn->right,i-lsize,delta);
	}
   }

	// access B[i], assumes i is right

static uint access (hybridBV B, uint64_t i, int64_t *delta, uint64_t n)

   { uint64_t lsize;
     if (B->type == tDynamic) 
        { B->bv.dyn->accesses++;
	  if (mustFlatten(B,n))
 	     flatten(B,delta); 
          else 
	     { lsize = hybridLength(B->bv.dyn->left);
               if (i < lsize) return access(B->bv.dyn->left,i,delta,n);
               else return access(B->bv.dyn->right,i-lsize,delta,n);
	     }
        }
     if (B->type == tLeaf) return leafAccess(B->bv.leaf,i);
     return staticAccess(B->bv.stat,i);
   }

uint hybridAccess (hybridBV B, uint64_t i)

   { int64_t delta = 0;
     uint64_t n = 0;
     if (B->type == tDynamic) n = B->bv.dyn->size;
     uint answ = access(B,i,&delta,n);
     if (delta) recompute(B,i,delta);
     return answ;
   }

        // read bits [i..i+l-1], onto D[j...]

static void sread (hybridBV B, uint64_t i, uint64_t l, uint64_t *D, uint64_t j,
		   uint *recomp, uint64_t n)

   { uint64_t lsize;
     int64_t delta;
     if (B->type == tDynamic)
        { B->bv.dyn->accesses++;
	  if (mustFlatten(B,n)) {
	     delta = 0;
	     flatten(B,&delta); 
	     if (delta) *recomp = 1;
	     }
	  else {
	    lsize = hybridLength(B->bv.dyn->left);
     	    if (i+l < lsize) sread(B->bv.dyn->left,i,l,D,j,recomp,n);
	    else if (i>=lsize) sread(B->bv.dyn->right,i-lsize,l,D,j,recomp,n);
	    else { sread(B->bv.dyn->left,i,lsize-i,D,j,recomp,n);
		   sread(B->bv.dyn->right,0,l-(lsize-i),D,j+(lsize-i),recomp,n);
		 }
	    return;
	    }
        }
     if (B->type == tLeaf) { leafRead(B->bv.leaf,i,l,D,j); return; }
     staticRead(B->bv.stat,i,l,D,j);
   }

void hybridRead (hybridBV B, uint64_t i, uint64_t l, uint64_t *D, uint64_t j)

   { uint recomp = 0;
     uint64_t n = 0;
     if (B->type == tDynamic) n = B->bv.dyn->size;
     sread(B,i,l,D,j,&recomp,n);
     if (recomp) rrecompute(B,i,l);
   }

	// computes rank_1(B,i), zero-based, assumes i is right

static uint64_t rank (hybridBV B, uint64_t i, int64_t *delta, uint64_t n)

   { uint64_t lsize;
     if (B->type == tDynamic)
        { B->bv.dyn->accesses++;
	  if (mustFlatten(B,n))
	     flatten(B,delta); 
          else { 
	     lsize = hybridLength(B->bv.dyn->left);
             if (i < lsize) return rank(B->bv.dyn->left,i,delta,n);
             else return hybridOnes(B->bv.dyn->left) + 
			 rank(B->bv.dyn->right,i-lsize,delta,n);
	     }
	}
     if (B->type == tLeaf) return leafRank(B->bv.leaf,i);
     return staticRank(B->bv.stat,i);
   }

uint64_t hybridRank (hybridBV B, uint64_t i)

   { int64_t delta = 0;
     uint64_t n = 0;
     if (B->type == tDynamic) n = B->bv.dyn->size;
     uint64_t answ = rank(B,i,&delta,n);
     if (delta) recompute(B,i,delta);
     return answ;
   }

	// computes rank_0(B,i), zero-based, assumes i is right

uint64_t hybridRank0 (hybridBV B, uint64_t i)

   { return i + 1 - hybridRank(B,i);
   }

        // computes select_1(B,j), zero-based, assumes j is right

static uint64_t select1 (hybridBV B, uint64_t j, int64_t *delta, uint64_t n)

   { uint64_t lones;
     if (B->type == tDynamic)
        { B->bv.dyn->accesses++;
	  if (mustFlatten(B,n))
	     flatten(B,delta); 
          else { 
             lones = hybridOnes(B->bv.dyn->left);
             if (j <= lones) return select1(B->bv.dyn->left,j,delta,n);
	     return hybridLength(B->bv.dyn->left) 
		    + select1(B->bv.dyn->right,j-lones,delta,n);
	     }
	}
     if (B->type == tLeaf) return leafSelect(B->bv.leaf,j);
     return staticSelect(B->bv.stat,j);
   }

uint64_t hybridSelect (hybridBV B, uint64_t j)

   { int64_t delta = 0;
     uint64_t n = 0;
     if (B->type == tDynamic) n = B->bv.dyn->size;
     uint64_t answ = select1(B,j,&delta,n);
     if (delta) recompute(B,answ,delta);
     return answ;
   }

        // computes select_0(B,j), zero-based, assumes j is right

static uint64_t select0 (hybridBV B, uint64_t j, int64_t *delta, uint64_t n)

   { uint64_t lzeros;
     if (B->type == tDynamic)
        { B->bv.dyn->accesses++;
	  if (mustFlatten(B,n))
	     flatten(B,delta); 
          else { 
             lzeros = hybridLength(B->bv.dyn->left)-hybridOnes(B->bv.dyn->left);
             if (j <= lzeros) return select0(B->bv.dyn->left,j,delta,n);
	     return hybridLength(B->bv.dyn->left) 
		    + select0(B->bv.dyn->right,j-lzeros,delta,n);
	     }
	}
     if (B->type == tLeaf) return leafSelect0(B->bv.leaf,j);
     return staticSelect0(B->bv.stat,j);
   }

uint64_t hybridSelect0 (hybridBV B, uint64_t j)

   { int64_t delta = 0;
     uint64_t n = 0;
     if (B->type == tDynamic) n = B->bv.dyn->size;
     uint64_t answ = select0(B,j,&delta,n);
     if (delta) recompute(B,answ,delta);
     return answ;
   }

        // computes next_1(B,i), zero-based and including i
        // returns -1 if no answer

static int64_t next1 (hybridBV B, uint64_t i, int64_t *delta, uint64_t n)

   { uint64_t lsize;
     int64_t next;
     if (B->type == tDynamic)
        { if (hybridOnes(B) == 0) return -1; // not considered an access!
	  B->bv.dyn->accesses++;
	  if (mustFlatten(B,n))
	     flatten(B,delta); 
          else { 
	     lsize = hybridLength(B->bv.dyn->left);
             if (i < lsize) 
		{ next = next1(B->bv.dyn->left,i,delta,n);
		  if (next != -1) return next;
		  i = lsize;
		}
	     next = next1(B->bv.dyn->right,i-lsize,delta,n);
	     if (next == -1) return -1;
             return lsize + next;
	     }
	}
     if (B->type == tLeaf) return leafNext(B->bv.leaf,i);
     return staticNext(B->bv.stat,i);
   }

int64_t hybridNext (hybridBV B, uint64_t i)

   { int64_t delta = 0;
     uint64_t n = 0;
     if (B->type == tDynamic) n = B->bv.dyn->size;
     int64_t answ = next1(B,i,&delta,n);
     if (delta) recompute(B,i,delta);
     return answ;
   }

        // computes next_0(B,i), zero-based and including i
        // returns -1 if no answer

static int64_t next0 (hybridBV B, uint64_t i, int64_t *delta, uint64_t n)

   { uint64_t lsize;
     int64_t next;
     if (B->type == tDynamic)
        { if (hybridOnes(B) == hybridLength(B)) return -1; // not an access
	  B->bv.dyn->accesses++;
	  if (mustFlatten(B,n))
	     flatten(B,delta); 
          else { 
	     lsize = hybridLength(B->bv.dyn->left);
             if (i < lsize) 
		{ next = next0(B->bv.dyn->left,i,delta,n);
		  if (next != -1) return next;
		  i = lsize;
		}
	     next = next0(B->bv.dyn->right,i-lsize,delta,n);
	     if (next == -1) return -1;
             return lsize + next;
	     }
	}
     if (B->type == tLeaf) return leafNext0(B->bv.leaf,i);
     return staticNext0(B->bv.stat,i);
   }

int64_t hybridNext0 (hybridBV B, uint64_t i)

   { int64_t delta = 0;
     uint64_t n = 0;
     if (B->type == tDynamic) n = B->bv.dyn->size;
     int64_t answ = next0(B,i,&delta,n);
     if (delta) recompute(B,i,delta);
     return answ;
   }

