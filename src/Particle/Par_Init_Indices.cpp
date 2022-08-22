#include "GAMER.h"

#ifdef PARTICLE

// this function pointer must be set by a test problem initializer
void (*Par_Init_Indices_Ptr)() = NULL;


void Par_Init_IndicesDefault()
{
   
   long myIndexOffset;

// Assign the particle indices on this processor
   for ( long p=0; p<amr->Par->NPar_AcPlusInac; p++ )
      amr->Par->Indx[p] = (real)p;

// Find the index offset for this processor
   Par_FindIndexOffset( amr->Par->NPar_AcPlusInac, &myIndexOffset );

// Add the offset to the particle indices to ensure global unique
// indices
   for ( long p=0; p<amr->Par->NPar_AcPlusInac; p++ )
      amr->Par->Indx[p] += (real)myIndexOffset;


} // FUNCTION : Par_Init_IndicesDefault

#endif // #ifdef PARTICLE