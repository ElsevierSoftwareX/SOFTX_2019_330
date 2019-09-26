#include <stdint.h>

union count {
 uint64_t count;
 uint64_t *p;
};

struct ht3{
 struct ht3e *contents;
 struct ht3e **map;
 int numElements;
};

struct ht3e{
 int abc[3];
 struct ht3e *nxt; //Next element of the same hash
 uint64_t count;
 int descCountOpPointerMask;
 union count desc[3];
};

// The idea goes like this
// 
// hte3~0 (a,1,!) count=#a1! desc=[ #a!, #2!, #! ]
// hte3~1 (a,2,!) count=#a2! desc=[ -> ~0.desc[0], #2!, -> ~0.desc[2] ]
// hte3~2 (b,1,@) count=#b1@ desc=[ #b@, #1@, #@ ]
// hte3~3 (a,1,@) count=#a1@ desc=[ #a@, -> ~2.desc[1], -> ~2.desc[2] ]
// ...
// 
// ht3 stores hash(X,Y,Z)->hte; many hashes can point to the same hte, it behaves like an index.
// Again, hash is (X+1)*(nZ+1)*(nY+1)+(Y+1)*(nZ+1)+(Z+1) when (nX+1)*(nZ+1)*(nY+1)+(nZ+1)*(nY+1)+(nZ+1)<=N (or N some hash index amplifier)
// or some true hash function of X,Y,Z otherwise.
// 
// Finding works like this (note that X,..,Z==-1 is valid and means who cares)
// 1. hash(X,Y,Z)->hte; if this hte.abc!=(X,Y,Z), go to next [-1==everything]. If next is NULL, say NOTFOUND
// 
// Adding works like this. It only works for X,Y,Z!=-1.
// 1. Find X,Y,Z. If found, increment count and each element in desc (if pointer, de-refence and then increment)
// 2. If not, add a new ht3e, linking last ht3e's next to it. 
// 3. Init count to 0.
// 4. For each desc, look for X,Y,Z with -1 in the element this desc ignores.
// 5. If found, investigate its corresponding dist copy its pointer or make a pointer to its count (IMO copying pointer will never occur)
// 6. If not found, init your desc as count=0, and add itself in a slot of this hash, by going on the ongoing next chain and replace final NULL with itself; still, do nothing if found itself on the chain.
