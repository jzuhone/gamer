#include "GAMER.h"
#include "BlackHole.h"



std::vector<BlackHole_t> CM_BH;



#if ( MODEL == HYDRO  &&  defined GRAVITY  &&  defined MASSIVE_PARTICLES )

//-------------------------------------------------------------------------------------------------------
// Function    :  BlackHole_t::UpdateAccretionRate
// Description :  Calculate and update the black hole accretion rate (Mdot_tot/Mdot_hot/Mdot_cold)
//
// Note        :  1. Uses the accretion-region diagnostics (RAcc_*) and the current mass already
//                   stored in this object
//
// Parameter   :  Mode  : Accretion mode (1: hot: 2: cold; 3: hot+cold)
//                r_acc : Accretion radius of the black hole
//
// Return      :  Mdot_tot, Mdot_hot, Mdot_cold (member variables)
//-------------------------------------------------------------------------------------------------------
void BlackHole_t::UpdateAccretionRate( const int Mode, const double r_acc )
{

   double acc_cold = 0.0, acc_hot = 0.0;

// hot accretion rate
   if (  ( Mode == 1 || Mode == 3 )  &&  RAcc_GasDens > 0.0  )
   {
      acc_hot = 4.0 * M_PI * SQR(NEWTON_G) * SQR(Mass) * RAcc_GasDens /
                pow( SQR(RAcc_SoundSpeed) + SQR(RAcc_RelativeVel), 1.5 );
   }

// cold accretion rate
   if (  ( Mode == 2 || Mode == 3 )  &&  (RAcc_GasMass+RAcc_ParMass) > 0.0  )
   {
      const double t_ff = sqrt( 2.0*CUBE(r_acc) / NEWTON_G / (RAcc_GasMass+RAcc_ParMass) );
      acc_cold = RAcc_ColdGasMass / t_ff;
   }

   Mdot_tot  = acc_hot + acc_cold;
   Mdot_hot  = acc_hot;
   Mdot_cold = acc_cold;

} // FUNCTION : BlackHole_t::UpdateAccretionRate



//-------------------------------------------------------------------------------------------------------
// Function    :  BlackHole_t::ResetSubstepCounters
// Description :  Reset the per-substep counters/quantities recomputed every call to
//                Flu_ResetByUser_API_ClusterMerger()
//-------------------------------------------------------------------------------------------------------
void BlackHole_t::ResetSubstepCounters( const double Jet_HalfHeight )
{

   Bondi_SinkNCell = 0;
   Jet_WaveK       = 0.5 * M_PI / Jet_HalfHeight;

} // FUNCTION : BlackHole_t::ResetSubstepCounters



//-------------------------------------------------------------------------------------------------------
// Function    :  BlackHole_t::ResetRecordDiagnostics
// Description :  Zero the cumulative sink/injection diagnostics after they have been written out by
//                Aux_Record_ClusterMerger()
//-------------------------------------------------------------------------------------------------------
void BlackHole_t::ResetRecordDiagnostics()
{

   Bondi_SinkMass    = 0.0;
   Bondi_SinkMomX    = 0.0;
   Bondi_SinkMomY    = 0.0;
   Bondi_SinkMomZ    = 0.0;
   Bondi_SinkMomXAbs = 0.0;
   Bondi_SinkMomYAbs = 0.0;
   Bondi_SinkMomZAbs = 0.0;
   Bondi_SinkE       = 0.0;
   Bondi_SinkEk      = 0.0;
   Bondi_SinkEt      = 0.0;
   E_inj_exp         = 0.0;
   M_inj_exp         = 0.0;

} // FUNCTION : BlackHole_t::ResetRecordDiagnostics



//-------------------------------------------------------------------------------------------------------
// Function    :  BlackHole_t::MergeFrom
// Description :  Absorb another black hole's mass into this one (used when two black holes merge)
//
// Parameter   :  Other : The black hole being merged into this one; its mass is zeroed out
//-------------------------------------------------------------------------------------------------------
void BlackHole_t::MergeFrom( BlackHole_t &Other )
{

   Mass       += Other.Mass;
   Other.Mass  = 0.0;

} // FUNCTION : BlackHole_t::MergeFrom

#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY  &&  defined MASSIVE_PARTICLES )
