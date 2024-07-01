
#ifndef DEFS
#define DEFS

#define FT_SIZE 13
#define FT_SIZE_WEIGHTS 15

#define GRAD_ELEMENT_SIZE 8

#define MAX(a,b) a>b?a:b

// current code assumes that we have 16 sketches - will not work otherwise
#define UnivMon_CS_LVLS 16

// what hash to use
//#define USE_XXHASH
//#define USE_BOBHASH
//#define USE_ONE_XXHASH
//#define USE_TWO_XXHASHES
#define USE_128_BIT_XXHASH


#endif // !DEFS


