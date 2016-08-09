#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_PassParticle2Son_AllPatch
// Description :  Pass particles from father to sons for all patches at the target level
//
// Note        :  1. It simply invokes Par_PassParticle2Son for all patches at the target level
//                   (and with sons at home)
//                   --> For patches with sons living abroad, this function will first collect
//                       particles from other ranks to the father-buffer patches in this rank
//                       by calling Par_LB_ExchangeParticleBetweenPatch, and then pass
//                       these particles to the real son patches in the same rank by calling
//                       Par_PassParticle2Son again
//                2. It is invoked in EvolveLevel after the velocity correction in KDK
//
// Parameter   :  FaLv  : Father's refinement level
//-------------------------------------------------------------------------------------------------------
void Par_PassParticle2Son_AllPatch( const int FaLv )
{

// nothing to do if there is no patch at FaLv+1
   if ( FaLv == TOP_LEVEL  ||  NPatchTotal[FaLv+1] == 0 )   return;


// pass particles from fathers to sons if they are in the same rank
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
   {
//    patches with sons in other ranks (those with SonPID<-1) will be dealt with later
      if ( amr->patch[0][FaLv][FaPID]->son >= 0 )  Par_PassParticle2Son( FaLv, FaPID );
   }


// deal with the case when fathers and sons are NOT in the same rank
#  ifdef LOAD_BALANCE
   int FaBufPID;

// collect particles from other ranks to the father-buffer patches in this rank
   Par_LB_ExchangeParticleBetweenPatch(
      FaLv,
      amr->Par->F2S_Send_NPatchTotal[FaLv], amr->Par->F2S_Send_PIDList[FaLv], amr->Par->F2S_Send_NPatchEachRank[FaLv],
      amr->Par->F2S_Recv_NPatchTotal[FaLv], amr->Par->F2S_Recv_PIDList[FaLv], amr->Par->F2S_Recv_NPatchEachRank[FaLv] );

// pass particles from father-buffer patches to their real son patches in the same rank
   for (int t=0; t<amr->Par->F2S_Recv_NPatchTotal[FaLv]; t++)
   {
      FaBufPID = amr->Par->F2S_Recv_PIDList[FaLv][t];

#     ifdef DEBUG_PARTICLE
      if ( FaBufPID < amr->NPatchComma[FaLv][1] )
         Aux_Error( ERROR_INFO, "This is NOT a father-buffer patch (FaLv %d, FaPID %d, NReal %d) !!\n",
                    FaLv, FaBufPID, amr->NPatchComma[FaLv][1] );

      if ( amr->patch[0][FaLv][FaBufPID]->son < 0 )
         Aux_Error( ERROR_INFO, "Father-buffer patch has no son at home (FaLv %d, FaPID %d, SonPID %d) !!\n",
                    FaLv, FaBufPID, amr->patch[0][FaLv][FaBufPID]->son );
#     endif

      Par_PassParticle2Son( FaLv, FaBufPID );
   }
#  endif // #ifdef LOAD_BALANCE

} // FUNCTION : Par_PassParticle2Son_AllPatch



#endif // #ifdef PARTICLE