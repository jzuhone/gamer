#include "GAMER.h"

#ifdef PARTICLE

void Par_FindIndexOffset( const long newCount, long *indexOffset )
{

    long startIndexNumber = amr->Par->StartIndexNumber;
    long *localNumParticles = new long [MPI_NRank];

    MPI_Allgather( &newCount, 1, MPI_LONG, localNumParticles, 
                   1, MPI_LONG, MPI_COMM_WORLD );

    long myIndexStart = 0;
    for ( int i = 0; i < MPI_NRank; i++ )
        myIndexStart += localNumParticles[i];
  
    *indexOffset = myIndexStart + startIndexNumber;

    delete [] localNumParticles;

//  Update value of startIndexNumber

    int maxIndexOwner = MPI_NRank-1; 

    if ( MPI_Rank == maxIndexOwner ) 
        startIndexNumber = *indexOffset + newCount;

    MPI_Bcast( startIndexNumber, 1, MPI_LONG;
               maxIndexOwner, MPI_COMM_WORLD );
    
    amr->Par->startIndexNumber = startIndexNumber;

    return;

} // FUNCTION : Par_FindIndexOffset

#endif // #ifdef PARTICLE