#include "GAMER.h"

// problem-specific global variables
// =======================================================================================
extern int      Merger_Coll_NumBHs;
extern double   R_acc;                // the radius to compute the accretion rate
extern double (*CM_ClusterCen)[3];
extern FieldIdx_t RefineFieldIdx;
// =======================================================================================



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_ClusterMerger
// Description :  Flag cells for refinement for the galaxy cluster merger simulation
//
// Note        :  1. Linked to the function pointer "Flag_User_Ptr" by "Init_TestProb_ClusterMerger()"
//                   to replace "Flag_User()"
//                2. Please turn on the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k       : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv          : Refinement level of the target patch
//                PID         : ID of the target patch
//                Threshold   : Useless here
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_ClusterMerger( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   const double dh     = amr->dh[lv];
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

   bool Flag1 = false;

// flag cells within the target radius, and if the radius is not resolved with a specific number (Threshold[0]) of cells
   for (int c=0; c<Merger_Coll_NumBHs; c++)
   {
      if ( DIST_SQR_3D( Pos, CM_ClusterCen[c] ) <= SQR(25*R_acc)  &&  R_acc/dh <= Threshold[0] )
      {
         Flag1 = true;
         break;
      } // if ( R_SQR <= SQR(25*R_acc)  &&  R_acc/dh <= Threshold[0] )
   } // for (int c=0; c<Merger_Coll_NumBHs; c++)

// flag if the mass fraction of the scalar exceeds the given threshold

   const real (*Rho )[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS];           // density
   const real (*Scal)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[RefineFieldIdx]; // passive scalar

   const real Frac = Scal[k][j][i] / Rho[k][j][i];
   bool Flag2 = Frac > Threshold[1];

   bool Flag = Flag1 || Flag2;

   return Flag;

} // FUNCTION : Flag_ClusterMerger
