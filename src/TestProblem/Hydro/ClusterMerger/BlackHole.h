#ifndef __CLUSTERMERGER_BLACKHOLE_H__
#define __CLUSTERMERGER_BLACKHOLE_H__



#include <vector>



//-------------------------------------------------------------------------------------------------------
// Structure   :  BlackHole_t
// Description :  Data structure storing the full state of a single black hole in the ClusterMerger
//                test problem (position, velocity, mass, accretion, jet feedback, and diagnostics)
//
// Note        :  1. One instance per black hole; all black holes are stored in the global vector CM_BH
//                   --> CM_BH is sized to Merger_Coll_NumHalos and never resized/erased; only indices
//                       [0, Merger_Coll_NumBHs) are considered "active" (a merger simply decrements
//                       Merger_Coll_NumBHs and leaves the merged-away slot stale, as before)
//                2. All members zero-initialize on construction
//-------------------------------------------------------------------------------------------------------
struct BlackHole_t
{

// position (also serves as the cluster center) and velocity
   double Pos[3] = { 0.0, 0.0, 0.0 };
   double Vel[3] = { 0.0, 0.0, 0.0 };
   double Mass   = 0.0;

// accretion rates
   double Mdot_tot = 0.0, Mdot_hot = 0.0, Mdot_cold = 0.0;

// jet geometry (static config loaded from Input__TestProb) and direction/feedback rates
   double Jet_HalfHeight = 0.0, Jet_Radius = 0.0;
   double Jet_Vec[3] = { 0.0, 0.0, 0.0 };
   double Jet_Mdot = 0.0, Jet_Pdot = 0.0, Jet_Edot = 0.0;

// transient per-substep injection quantities, recomputed every call to Flu_ResetByUser_API_ClusterMerger()
   double Jet_WaveK = 0.0, V_cyl = 0.0, M_inj = 0.0, P_inj = 0.0, E_inj = 0.0, normalize_const = 0.0;

// accretion-region diagnostics
   double RAcc_GasVel[3] = { 0.0, 0.0, 0.0 };
   double RAcc_SoundSpeed = 0.0, RAcc_GasDens = 0.0, RAcc_RelativeVel = 0.0;
   double RAcc_ColdGasMass = 0.0, RAcc_GasMass = 0.0, RAcc_ParMass = 0.0;
   int    Cluster_NPar_close = 0;

// sink-region diagnostics, accumulated across sub-steps and drained by Aux_Record_ClusterMerger()
   double Bondi_SinkMass = 0.0, Bondi_SinkMomX = 0.0, Bondi_SinkMomY = 0.0, Bondi_SinkMomZ = 0.0;
   double Bondi_SinkMomXAbs = 0.0, Bondi_SinkMomYAbs = 0.0, Bondi_SinkMomZAbs = 0.0;
   double Bondi_SinkE = 0.0, Bondi_SinkEk = 0.0, Bondi_SinkEt = 0.0;
   int    Bondi_SinkNCell = 0;

// cumulative expected injection, also drained by Aux_Record_ClusterMerger()
   double E_inj_exp = 0.0, M_inj_exp = 0.0;


   void UpdateAccretionRate( const int Mode, const double r_acc );
   void ResetSubstepCounters();
   void ResetRecordDiagnostics();
   void MergeFrom( BlackHole_t &Other );

}; // struct BlackHole_t



extern std::vector<BlackHole_t> CM_BH;



#endif // #ifndef __CLUSTERMERGER_BLACKHOLE_H__
