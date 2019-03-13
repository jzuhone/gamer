#include "GAMER.h"
#include <typeinfo>

#ifdef VISCOSITY

//-------------------------------------------------------------------------------------------------------
// Function    :  Viscosity_Init
// Description :  Initialize viscosity
//
// Note        :  1. Must be called AFTER Init_Load_Parameter() and Init_Unit()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Viscosity_Init()
{

  if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

  if ( VISCOSITY_TYPE == SPITZER_VISCOSITY ) {
    real lnLambda = 30; // Coulomb logarithm for Spitzer viscosity
    // This prefactor is in CGS and we must convert it to code units
    // To avoid precision errors, we do this one step at a time
    // Still need to correct for differences in mean molecular weight
    FreqPrefactor = 2.68273501e34; // in cm**6/(g*s**4)
    FreqPrefactor *= UNIT_D;
    FreqPrefactor *= UNIT_D;
    FreqPrefactor /= UNIT_M;
    FreqPrefactor *= UNIT_T;
    FreqPrefactor *= UNIT_T;
    FreqPrefactor *= UNIT_T;
    FreqPrefactor *= UNIT_T;
    FreqPrefactor *= lnLambda;
  }

} // FUNCTION : Viscosity_Init

#endif // #ifdef VISCOSITY
