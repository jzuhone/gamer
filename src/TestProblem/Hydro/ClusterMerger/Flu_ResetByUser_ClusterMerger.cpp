#include "GAMER.h"
#include "BlackHole.h"



#if ( MODEL == HYDRO  &&  defined GRAVITY  &&  defined MASSIVE_PARTICLES )
// problem-specific global variables
// =======================================================================================
// (1) read from Input__TestProb
extern int        Merger_Coll_NumHalos;

extern bool       AGN_feedback;
extern int        Accretion_Mode, JetDirection_case;
extern double     eta, eps_f, eps_m, R_acc, R_dep;
extern double     Jet_HalfHeight, Jet_Radius;
extern bool       fixBH, AdjustBHPos, AdjustBHVel;
extern double     AdjustPeriod;
// ---------------------------------------------------------------------------------------
// (2) read from files
extern int        JetDirection_NBin;
extern double    *CM_Jet_Time_table;
extern double   **CM_Jet_Theta_table;
extern double   **CM_Jet_Phi_table;
// ---------------------------------------------------------------------------------------
// (3) other variables
extern int        AdjustCount;
extern int        Merger_Coll_NumBHs;

extern FieldIdx_t Idx_ParHalo;

extern long_par  *CM_ClusterIdx_Cur;

static bool       if_overlap = false;  // variable to keep track if the jet injection regions are overlapping
// =======================================================================================


// problem-specific function prototypes
// =======================================================================================
// (1) external functions
#ifdef MHD
extern double (*MHD_ResetByUser_VecPot_Ptr)( const double x, const double y, const double z, const double Time,
                                             const double dt, const int lv, const char Component, double AuxArray[] );
extern double (*MHD_ResetByUser_BField_Ptr)( const double x, const double y, const double z, const double Time,
                                             const double dt, const int lv, const char Component, double AuxArray[],
                                             const double B_in, const bool UseVecPot, const real *Ax, const real *Ay,
                                             const real *Az, const int i, const int j, const int k );
#endif
// ---------------------------------------------------------------------------------------
// (2) internal functions
static void GetClusterCenter( const int lv, const bool AdjustPos, const bool AdjustVel, std::vector<BlackHole_t> &BH );
static void SetJetDirection( const double TimeNew, const int lv, const int FluSg );
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Func_ClusterMerger
// Description :  Function to reset the fluid field in the cluster merger problem
//
// Note        :  1. Invoked by Flu_ResetByUser_API_ClusterMerger() and Hydro_Init_ByFunction_AssignData() using the
//                   function pointer "Flu_ResetByUser_Func_Ptr"
//                   --> This function pointer is reset by Init_TestProb_Hydro_ClusterMerger()
//                   --> Hydro_Init_ByFunction_AssignData(): constructing initial condition
//                       Flu_ResetByUser_API_ClusterMerger(): after each update
//                2. Input fluid[] stores the original values
//                3. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically
//                4. Enabled by the runtime option "OPT__RESET_FLUID"
//
// Parameter   :  fluid    : Fluid array storing both the input (original) and reset values
//                           --> Including both active and passive variables
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                dt       : Time interval to advance solution
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  0 : This cell has not been reset
//                1 : This cell has been reset by cluster 1
//                2 : This cell has been reset by cluster 2
//-------------------------------------------------------------------------------------------------------
int Flu_ResetByUser_Func_ClusterMerger( real fluid[], const double Emag, const double x, const double y, const double z,
                                        const double Time, const double dt, const int lv, double AuxArray[] )
{

   const double Pos[3] = { x, y, z };

// (1) SMBH Accretion
   const double V_dep = 4.0 / 3.0 * M_PI * CUBE( R_dep ); // the volume to remove gas
   double D_dep[Merger_Coll_NumBHs]; // the density need to be removed

   for (int c=0; c<Merger_Coll_NumBHs; c++)
   {
      D_dep[c] = CM_BH[c].Mdot_tot*dt/V_dep;
      double r2 = 0.0;

      for (int d=0; d<3; d++)   r2 += SQR( Pos[d] - CM_BH[c].Pos[d] );

      if ( r2 <= SQR(R_dep) )
      {
         const real dens_old = fluid[DENS];
         fluid[DENS] -= (real)D_dep[c];

         if ( fluid[DENS] < MIN_DENS )  fluid[DENS] = MIN_DENS;

         const real DepRatio = fluid[DENS] / dens_old;
         fluid[MOMX] *= DepRatio;
         fluid[MOMY] *= DepRatio;
         fluid[MOMZ] *= DepRatio;
         fluid[ENGY] *= DepRatio;
      } // if ( r2[c] <= SQR(R_dep) )
   } // for (int c=0; c<Merger_Coll_NumBHs; c++)


// (2) Jet Feedback
   double Jet_dr, Jet_dh;
   double Dis_c2m, Vec_c2m[3];
   real   EngySin;
   int    status_overlap, n_jet;
   int    which_cluster = 0;

   if ( if_overlap )
   {
      if ( Merger_Coll_NumBHs == 1 )   Aux_Error( ERROR_INFO, "Merger_Coll_NumBHs = 1 but if_overlap = true!\n" );

      if ( CM_BH[0].Jet_Edot >= CM_BH[1].Jet_Edot )   status_overlap = 0;   // only inject cluster 1
      else                                            status_overlap = 1;   // only inject cluster 2

      n_jet = Merger_Coll_NumBHs - 1;
   }

   else // if ( if_overlap )
   {
      status_overlap = 0;
      n_jet = Merger_Coll_NumBHs;
   } // if ( if_overlap ) ... else ...

   for (int c=status_overlap; c<(n_jet+status_overlap); c++)
   {
//    distance: jet center to mesh
      for (int d=0; d<3; d++)   Vec_c2m[d] = Pos[d] - CM_BH[c].Pos[d];

      Dis_c2m = sqrt( SQR(Vec_c2m[0]) + SQR(Vec_c2m[1]) + SQR(Vec_c2m[2]) );

      Jet_dh = fabs( CM_BH[c].Jet_Vec[0]*Vec_c2m[0] + CM_BH[c].Jet_Vec[1]*Vec_c2m[1] + CM_BH[c].Jet_Vec[2]*Vec_c2m[2] );
      Jet_dr = sqrt( SQR(Dis_c2m) - SQR(Jet_dh) );

      if ( Jet_dh <= Jet_HalfHeight  &&  Jet_dr <= Jet_Radius )
      {
         which_cluster += c+1;

         const bool Check_MinEint_No = false;

//       record the old fluid variables
         const real dens_old = fluid[DENS];
         const real eint_old = Hydro_Con2Eint( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ],
                                               fluid[ENGY], Check_MinEint_No, NULL_REAL, PassiveFloorMask, Emag,
                                               NULL, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                               h_EoS_Table );

//       accrete mass
         fluid[DENS] += CM_BH[c].M_inj;

//       transfer into BH frame
         fluid[MOMX] -= CM_BH[c].Vel[0]*dens_old;
         fluid[MOMY] -= CM_BH[c].Vel[1]*dens_old;
         fluid[MOMZ] -= CM_BH[c].Vel[2]*dens_old;

//       use a sine function to make the velocity smooth within the jet from +CM_BH[c].Jet_Vec to -CM_BH[c].Jet_Vec
         EngySin = CM_BH[c].E_inj*CM_BH[c].normalize_const*sin( CM_BH[c].Jet_WaveK*Jet_dh );

//       the new momentum is calculated from the old density, new density, old momentum and injected energy
//       the momentum injection alters only the component along the jet direction
         const real P_old_sqr = SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]);
         const real Ekin_old  = ( dens_old == (real)0.0 ) ? (real)0.0 : (real)0.5*P_old_sqr/dens_old;
         const real P_new     = SQRT( (real)2.0*fluid[DENS]*(EngySin + Ekin_old) );
         const real JetSign   = SIGN( Vec_c2m[0]*CM_BH[c].Jet_Vec[0] + Vec_c2m[1]*CM_BH[c].Jet_Vec[1] + Vec_c2m[2]*CM_BH[c].Jet_Vec[2] );

         real P_old_perp[3], P_old_para, P_new_para, P_new_para_sqr;
         P_old_para     = fluid[MOMX]*CM_BH[c].Jet_Vec[0] + fluid[MOMY]*CM_BH[c].Jet_Vec[1] + fluid[MOMZ]*CM_BH[c].Jet_Vec[2];
         P_old_perp[0]  = fluid[MOMX] - P_old_para*CM_BH[c].Jet_Vec[0];
         P_old_perp[1]  = fluid[MOMY] - P_old_para*CM_BH[c].Jet_Vec[1];
         P_old_perp[2]  = fluid[MOMZ] - P_old_para*CM_BH[c].Jet_Vec[2];
         P_new_para_sqr = SQR(P_new) - SQR(P_old_perp[0]) - SQR(P_old_perp[1]) - SQR(P_old_perp[2]);
         P_new_para_sqr = FMAX( P_new_para_sqr, (real)0.0 );   // avoid negative squared values caused by rounding errors
         P_new_para     = SQRT( P_new_para_sqr );
         P_new_para    *= JetSign;

         fluid[MOMX] = P_new_para*CM_BH[c].Jet_Vec[0] + P_old_perp[0];
         fluid[MOMY] = P_new_para*CM_BH[c].Jet_Vec[1] + P_old_perp[1];
         fluid[MOMZ] = P_new_para*CM_BH[c].Jet_Vec[2] + P_old_perp[2];

//       transfer back into the rest frame
         fluid[MOMX] += CM_BH[c].Vel[0]*fluid[DENS];
         fluid[MOMY] += CM_BH[c].Vel[1]*fluid[DENS];
         fluid[MOMZ] += CM_BH[c].Vel[2]*fluid[DENS];

         fluid[ENGY]  = (real)0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS] +
                        eint_old + Emag;
      } // if ( Jet_dh <= Jet_HalfHeight  &&  Jet_dr <= Jet_Radius )
   } // for (int c=status_overlap; c<(n_jet+status_overlap); c++)

   if ( which_cluster >= 3 )   Aux_Error( ERROR_INFO, "which_cluster (%d) >= 3!\n", which_cluster );

   return which_cluster;

} // FUNCTION : Flu_ResetByUser_Func_ClusterMerger



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_API_ClusterMerger
// Description :  API for resetting the fluid array in the Bondi accretion problem
//
// Note        :  1. Enabled by the runtime option "OPT__RESET_FLUID"
//                2. Invoked using the function pointer "Flu_ResetByUser_API_Ptr"
//                   --> This function pointer is reset by Init_TestProb_Hydro_ClusterMerger()
//                3. Currently does not work with "OPT__OVERLAP_MPI"
//                4. Invoke Flu_ResetByUser_Func_ClusterMerger() directly
//                5. Currently, Jet_Pdot (the jet change in momentum), P_inj (the injected
//                   momentum) and eps_m (the fraction of total energy that goes into the thermal
//                   energy) do not participate in any calculations.
//
// Parameter   :  lv    : Target refinement level
//                FluSg : Target fluid sandglass
//                MagSg   : Target B field sandglass
//                TimeNew : Current physical time (system has been updated from TimeOld to TimeNew in EvolveLevel())
//                dt      : Time interval to advance solution (can be different from TimeNew-TimeOld in COMOVING)
//-------------------------------------------------------------------------------------------------------
void Flu_ResetByUser_API_ClusterMerger( const int lv, const int FluSg, const int MagSg, const double TimeNew, const double dt )
{

   const bool    CurrentMaxLv = (  NPatchTotal[lv] > 0  &&  ( lv == MAX_LEVEL || NPatchTotal[lv+1] == 0 )  );
   const double  dh           = amr->dh[lv];
   const double  dv           = CUBE(dh);

// (1) get the BH position and velocity and adjust them if needed
   bool AdjustPosNow = false;
   bool AdjustVelNow = false;

// only adjust the BHs on the current maximum level
   if ( CurrentMaxLv  &&  AdjustCount < int(TimeNew/AdjustPeriod) )
   {
      if ( AdjustBHPos == true )   AdjustPosNow = true;
      if ( AdjustBHVel == true )   AdjustVelNow = true;
      AdjustCount += 1;
   } // if ( CurrentMaxLv  &&  AdjustCount < int(TimeNew/AdjustPeriod))

   GetClusterCenter( lv, AdjustPosNow, AdjustVelNow, CM_BH );


// (2) decide whether to merge BHs
   if ( Merger_Coll_NumBHs == 2 )
   {
      const double AbsRelPos = DIST_3D_DBL( CM_BH[0].Pos, CM_BH[1].Pos );
      const double AbsRelVel = DIST_3D_DBL( CM_BH[0].Vel, CM_BH[1].Vel );
      const double soften    = amr->dh[MAX_LEVEL];
      double escape_vel;

      if ( AbsRelPos > soften ) {
         escape_vel = sqrt( 2.0 * NEWTON_G * (CM_BH[0].Mass+CM_BH[1].Mass) / AbsRelPos );
      }

      else {
         escape_vel = sqrt( 2.0 * NEWTON_G * (CM_BH[0].Mass+CM_BH[1].Mass) / soften );
      } // if ( AbsRelPos > soften ) ... else ...

//    merge the two BHs if they are located within R_acc, and the relative velocity is small enough
      if ( AbsRelPos < R_acc  &&  AbsRelVel < 3.0*escape_vel )
      {
         int merge_index = 0;    // record BH 1 merge BH 2 / BH 2 merge BH 1
         Merger_Coll_NumBHs -= 1;
         if ( CM_BH[0].Mass >= CM_BH[1].Mass )   merge_index = 1;
         else                                    merge_index = 2;
         CM_BH[0].MergeFrom( CM_BH[1] );

//       relabel the BH and DM particles being merged
//       NOTE: this assumes that BHs are only supported in a binary merger
         CM_ClusterIdx_Cur[1] = 0; // The BH from cluster 1 now belongs to cluster 0
         for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
         {
            if ( amr->Par->AttributeInt[Idx_ParHalo][p] == (long_par)1  &&
                 amr->Par->Mass[p] >= (real_par)0.0  &&
                 amr->Par->Type[p] == PTYPE_BLACK_HOLE )    amr->Par->Type[p] = PTYPE_DARK_MATTER;
         } // for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)

         if ( MPI_Rank == 0 )
         Aux_Message( stdout, "BHs Merge! TimeNew = %14.8e; merge_index = %d, "
                              "BHPos1 = %14.8e, %14.8e, %14.8e; BHPos2 = %14.8e, %14.8e, %14.8e; "
                              "BHVel1 = %14.8e, %14.8e, %14.8e; BHVel2 = %14.8e, %14.8e, %14.8e; "
                              "AbsRelPos = %14.8e, AbsRelVel = %14.8e, escape_vel = %14.8e.\n",
                      TimeNew, merge_index,
                      CM_BH[0].Pos[0], CM_BH[0].Pos[1], CM_BH[0].Pos[2], CM_BH[1].Pos[0], CM_BH[1].Pos[1], CM_BH[1].Pos[2],
                      CM_BH[0].Vel[0], CM_BH[0].Vel[1], CM_BH[0].Vel[2], CM_BH[1].Vel[0], CM_BH[1].Vel[1], CM_BH[1].Vel[2],
                      AbsRelPos, AbsRelVel, escape_vel );
      } // if ( AbsRelPos < R_acc  &&  AbsRelVel < 3*escape_vel )
   } // if ( Merger_Coll_NumBHs == 2 )


   if ( AGN_feedback )
   {
//    (3) set the jet direction vector
      SetJetDirection( TimeNew, lv, FluSg );


//    (4) calculate the accretion and feedback
//    reset to 0 since we only want to record the number of void cells **for one sub-step**
      for (int c=0; c<Merger_Coll_NumBHs; c++)   CM_BH[c].ResetSubstepCounters( Jet_HalfHeight );

//    variables for each rank
      int    num_hot    [Merger_Coll_NumBHs];     // the number of hot gas cells inside the accretion radius
      double gas_mass   [Merger_Coll_NumBHs];     // total gas mass inside the accretion radius
      double rho_hot    [Merger_Coll_NumBHs];     // the average density inside the accretion radius (hot gas)
      double mass_cold  [Merger_Coll_NumBHs];     // cold gas mass (T < 5e5 K) inside the accretion radius
      double Cs_hot     [Merger_Coll_NumBHs];     // the average sound speed inside the accretion radius (hot gas)
      double gas_mom_hot[Merger_Coll_NumBHs][3];  // average hot gas momentum
      double V_cyl_exact[Merger_Coll_NumBHs];     // exact volume of jet cylinder
      double normalize  [Merger_Coll_NumBHs];     // for computing the correct normalization constant
      bool   if_overlap_each_rank = false;        // check for jet injection site overlap on current MPI rank

      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         num_hot    [c]    = 0;
         gas_mass   [c]    = 0.0;
         rho_hot    [c]    = 0.0;
         mass_cold  [c]    = 0.0;
         Cs_hot     [c]    = 0.0;
         V_cyl_exact[c]    = 0.0;
         normalize  [c]    = 0.0;
         for (int d=0; d<3; d++)
         gas_mom_hot[c][d] = 0.0;
      }

//    variables for all ranks
      int    num_hot_sum[Merger_Coll_NumBHs];
      double rho_hot_sum[Merger_Coll_NumBHs], Cs_hot_sum[Merger_Coll_NumBHs], gas_mom_hot_sum[Merger_Coll_NumBHs][3];
      double gas_vel_hot_sum[Merger_Coll_NumBHs][3], V_cyl_exact_sum[Merger_Coll_NumBHs], normalize_sum[Merger_Coll_NumBHs];

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

         for (int k=0; k<PS1; k++)   { const double z = z0 + k*dh;
         for (int j=0; j<PS1; j++)   { const double y = y0 + j*dh;
         for (int i=0; i<PS1; i++)   { const double x = x0 + i*dh;

            const double Pos[3] = { x, y, z };
            real fluid_acc[NCOMP_TOTAL];

            for (int v=0; v<NCOMP_TOTAL; v++)   fluid_acc[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

#           ifdef MHD
            const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
#           else
            const real Emag = (real)0.0;
#           endif

            int status_overlap = 0;

            for (int c=0; c<Merger_Coll_NumBHs; c++)
            {
//             calculate the average density, sound speed and gas velocity inside accretion radius
               if (  DIST_SQR_3D( Pos, CM_BH[c].Pos ) <= SQR(R_acc)  )
               {
                  gas_mass[c] += fluid_acc[0]*dv;

#                 ifdef DUAL_ENERGY
                  const real Pres = Hydro_DensDual2Pres( fluid_acc[DENS], fluid_acc[DUAL], EoS_AuxArray_Flt[1],
                                                         false, NULL_REAL );
                  const real Eint = EoS_DensPres2Eint_CPUPtr( fluid_acc[DENS], Pres, NULL, EoS_AuxArray_Flt,
                                                              EoS_AuxArray_Int, h_EoS_Table );
                  const real Temp = EoS_DensEint2Temp_CPUPtr( fluid_acc[DENS], Eint, NULL, EoS_AuxArray_Flt,
                                                              EoS_AuxArray_Int, h_EoS_Table );
#                 else
                  const real Pres = Hydro_Con2Pres( fluid_acc[DENS], fluid_acc[MOMX], fluid_acc[MOMY], fluid_acc[MOMZ],
                                                    fluid_acc[ENGY], fluid_acc+NCOMP_FLUID, true, MIN_PRES, PassiveFloorMask, Emag,
                                                    EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr,
                                                    EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                    h_EoS_Table, NULL );
                  const real Temp = Hydro_Con2Temp( fluid_acc[DENS], fluid_acc[MOMX], fluid_acc[MOMY], fluid_acc[MOMZ],
                                                    fluid_acc[ENGY], fluid_acc+NCOMP_FLUID, true, MIN_TEMP, PassiveFloorMask, Emag,
                                                    EoS_DensEint2Temp_CPUPtr, EoS_GuessHTilde_CPUPtr,
                                                    EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                    h_EoS_Table );
#                 endif

                  const double cold_temp_thresh = 5e5;
                  if ( Temp <= cold_temp_thresh )   mass_cold[c] += fluid_acc[DENS]*dv;
                  else
                  {
                     rho_hot[c] += fluid_acc[DENS]*dv;
//                   NOTE: currently, the average sound speed is computed without applying any weighting.
                     Cs_hot[c] += sqrt( EoS_DensPres2CSqr_CPUPtr( fluid_acc[DENS], Pres, NULL, EoS_AuxArray_Flt,
                                        EoS_AuxArray_Int, h_EoS_Table ) );
                     for (int d=0; d<3; d++)   gas_mom_hot[c][d] += fluid_acc[d+MOMX]*dv;
                     num_hot[c] += 1;
                  } // if ( Temp <= cold_temp_thresh ) ... else ...
               } // if (  DIST_SQR_3D( Pos, CM_BH[c].Pos ) <= SQR(R_acc)  )

//             calculate the exact volume of jet cylinder and normalization
               if ( CurrentMaxLv )
               {
                  double Jet_dr, Jet_dh, Dis_c2m, Vec_c2m[3];

                  for (int d=0; d<3; d++)   Vec_c2m[d] = Pos[d] - CM_BH[c].Pos[d];
                  Dis_c2m = sqrt( SQR(Vec_c2m[0]) + SQR(Vec_c2m[1]) + SQR(Vec_c2m[2]) );
                  Jet_dh = fabs( CM_BH[c].Jet_Vec[0]*Vec_c2m[0] + CM_BH[c].Jet_Vec[1]*Vec_c2m[1] + CM_BH[c].Jet_Vec[2]*Vec_c2m[2] );
                  Jet_dr = sqrt( SQR(Dis_c2m) - SQR(Jet_dh) );

                  if ( Jet_dh <= Jet_HalfHeight  &&  Jet_dr <= Jet_Radius )
                  {
                     V_cyl_exact[c] += dv;
                     normalize[c]   += sin( CM_BH[c].Jet_WaveK*Jet_dh )*dv;
                     status_overlap ++;
                  } // if ( Jet_dh <= Jet_HalfHeight  &&  Jet_dr <= Jet_Radius )
               } // if ( CurrentMaxLv )
            } // for (int c=0; c<Merger_Coll_NumBHs; c++)

            if ( status_overlap >= 2 )   if_overlap_each_rank = true;

         }}} // for (int k=0; k<PS1; k++); for (int j=0; j<PS1; j++); for (int i=0; i<PS1; i++)
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         MPI_Allreduce( &num_hot    [c], &num_hot_sum        [c], 1, MPI_INT,    MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &gas_mass   [c], &CM_BH[c].RAcc_GasMass,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &rho_hot    [c], &rho_hot_sum        [c], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &mass_cold  [c], &CM_BH[c].RAcc_ColdGasMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &Cs_hot     [c], &Cs_hot_sum         [c], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( gas_mom_hot [c],  gas_mom_hot_sum    [c], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &V_cyl_exact[c], &V_cyl_exact_sum    [c], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         MPI_Allreduce( &normalize  [c], &normalize_sum      [c], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      } // for (int c=0; c<Merger_Coll_NumBHs; c++)
      MPI_Allreduce( &if_overlap_each_rank, &if_overlap, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD );

//    calculate the total DM mass inside the accretion radius
      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         if ( Accretion_Mode == 2  ||  Accretion_Mode == 3 )
         {
            double par_mass = 0.0;
            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            {
               const double *EdgeL        = amr->patch[0][lv][PID]->EdgeL;
               const double *EdgeR        = amr->patch[0][lv][PID]->EdgeR;
               const double  patch_pos[3] = { (EdgeL[0]+EdgeR[0])*0.5, (EdgeL[1]+EdgeR[1])*0.5, (EdgeL[2]+EdgeR[2])*0.5 };
               const double  patch_d      = sqrt( SQR(EdgeL[0]-EdgeR[0]) + SQR(EdgeL[1]-EdgeR[1]) + SQR(EdgeL[2]-EdgeR[2]) ) * 0.5;

               if (  DIST_SQR_3D( patch_pos, CM_BH[c].Pos ) <= SQR(2*R_acc+patch_d)  )
               for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               {
                  const long     ParID     = amr->patch[0][lv][PID]->ParList[p];
                  const real_par ParM      = amr->Par->Mass[ParID];
                  const real_par ParPos[3] = { amr->Par->PosX[ParID], amr->Par->PosY[ParID], amr->Par->PosZ[ParID] };

                  if (  DIST_SQR_3D( ParPos, CM_BH[c].Pos ) <= SQR(R_acc)  )    par_mass += ParM;
               } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
            } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

            MPI_Allreduce( &par_mass, &CM_BH[c].RAcc_ParMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         } // if ( Accretion_Mode == 2  ||  Accretion_Mode == 3 )
      } // for (int c=0; c<Merger_Coll_NumBHs; c++)

//    calculate the accretion rate
      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         if ( num_hot_sum[c] == 0 )
         {
            CM_BH[c].RAcc_GasDens     = 0.0;
            CM_BH[c].RAcc_SoundSpeed  = 0.0;
            CM_BH[c].RAcc_RelativeVel = 0.0;
            for (int d=0; d<3; d++)
            CM_BH[c].RAcc_GasVel[d]   = 0.0;
         }

         else
         {
            double RelVelSqr = 0.0;  // the relative velocity between BH and gas
            for (int d=0; d<3; d++)  gas_vel_hot_sum[c][d] = gas_mom_hot_sum[c][d] / rho_hot_sum[c];
//###ISSUE: Explore using V_cyl_exact_sum[c] as the volume
            rho_hot_sum[c] /= ( 4.0 / 3.0 * M_PI * CUBE(R_acc) );
            Cs_hot_sum[c]  /= (double)num_hot_sum[c];
            for (int d=0; d<3; d++)  RelVelSqr += SQR( CM_BH[c].Vel[d] - gas_vel_hot_sum[c][d] );

            CM_BH[c].RAcc_GasDens     = rho_hot_sum    [c];
            CM_BH[c].RAcc_SoundSpeed  = Cs_hot_sum     [c];
            CM_BH[c].RAcc_RelativeVel = sqrt( RelVelSqr );
            for (int d=0; d<3; d++)
            CM_BH[c].RAcc_GasVel[d]   = gas_vel_hot_sum[c][d];
         } // if ( num_hot_sum[c] == 0 ) ... else ...

         CM_BH[c].UpdateAccretionRate( Accretion_Mode, R_acc );

         if ( V_cyl_exact_sum[c] != 0 )   CM_BH[c].normalize_const = V_cyl_exact_sum[c] / normalize_sum[c];
         else                             CM_BH[c].normalize_const = 0.5 * M_PI;
      } // for (int c=0; c<Merger_Coll_NumBHs; c++)

//    update BH mass
      for (int c=0; c<Merger_Coll_NumBHs; c++)   if ( CurrentMaxLv )   CM_BH[c].Mass += CM_BH[c].Mdot_tot * dt;


//    (5) calculate the injection rate
      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         CM_BH[c].Jet_Mdot = eta * CM_BH[c].Mdot_tot;
         CM_BH[c].Jet_Pdot = sqrt( 2.0*eta*eps_f*(1.0-eps_m) ) * CM_BH[c].Mdot_tot * (Const_c/UNIT_V);
         CM_BH[c].Jet_Edot = eps_f * CM_BH[c].Mdot_tot * SQR( Const_c/UNIT_V );
         CM_BH[c].V_cyl    = M_PI * SQR(Jet_Radius) * 2.0 * Jet_HalfHeight;

//       calculate the density that needs to be injected
         if ( CurrentMaxLv  &&  V_cyl_exact_sum[c] != 0.0 )
         {
            CM_BH[c].M_inj = CM_BH[c].Jet_Mdot * dt / V_cyl_exact_sum[c];
            CM_BH[c].P_inj = CM_BH[c].Jet_Pdot * dt / V_cyl_exact_sum[c];
            CM_BH[c].E_inj = CM_BH[c].Jet_Edot * dt / V_cyl_exact_sum[c];
         }

         else
         {
            CM_BH[c].M_inj = CM_BH[c].Jet_Mdot * dt / CM_BH[c].V_cyl;
            CM_BH[c].P_inj = CM_BH[c].Jet_Pdot * dt / CM_BH[c].V_cyl;
            CM_BH[c].E_inj = CM_BH[c].Jet_Edot * dt / CM_BH[c].V_cyl;
         } // if ( CurrentMaxLv  &&  V_cyl_exact_sum[c] != 0.0 ) ... else ...
      } // for (int c=0; c<Merger_Coll_NumBHs; c++)

      if ( CurrentMaxLv )
      {
         for (int c=0; c<Merger_Coll_NumBHs; c++)   CM_BH[c].E_inj_exp += CM_BH[c].Jet_Edot*dt;
         for (int c=0; c<Merger_Coll_NumBHs; c++)   CM_BH[c].M_inj_exp += CM_BH[c].Jet_Mdot*dt;
      } // if ( CurrentMaxLv )


//    (6) perform injection
//    use the "static" schedule for reproducibility
#     pragma omp parallel for schedule( static )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         int Reset = 0;
         const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

         for (int k=0; k<PS1; k++)  {  const double z = z0 + k*dh;
         for (int j=0; j<PS1; j++)  {  const double y = y0 + j*dh;
         for (int i=0; i<PS1; i++)  {  const double x = x0 + i*dh;
            real fluid[NCOMP_TOTAL], fluid_old[NCOMP_TOTAL];

            for (int v=0; v<NCOMP_TOTAL; v++)
            {
               fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

//             backup the unmodified values since we want to record the amount of sunk variables removed at the maximum level
               fluid_old[v] = fluid[v];
            }

#           ifdef MHD
            const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
#           else
            const real Emag = (real)0.0;
#           endif

//          reset this cell
            if ( CurrentMaxLv )   Reset = Flu_ResetByUser_Func_ClusterMerger( fluid, Emag, x, y, z, TimeNew, dt, lv, NULL );

//          operations necessary only when this cell has been reset
            if ( Reset != 0 )
            {
#              if ( MODEL == HYDRO  &&  !defined SRHD )
//             abort the program instead of applying a density floor
               if ( fluid[DENS] < MIN_DENS )   Aux_Error( ERROR_INFO, "Fluid density has reached the floor!\n" );

//             apply an energy floor
               fluid[ENGY] = Hydro_CheckMinEintInEngy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                                       (real)MIN_EINT, PassiveFloorMask, Emag );

//             calculate the dual-energy variable
#              ifdef DUAL_ENERGY
               fluid[DUAL] = Hydro_Con2Dual( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], fluid+NCOMP_FLUID, Emag,
                                             EoS_DensEint2Pres_CPUPtr, EoS_CREint2CRPres_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                             PassiveFloorMask );
#              endif

//             floor and normalize passive scalars
#              if ( NCOMP_PASSIVE > 0 )
               for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)   fluid[v] = FMAX( fluid[v], TINY_NUMBER );

               if ( OPT__NORMALIZE_PASSIVE )
                  Hydro_NormalizePassive( fluid[DENS], fluid+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#              endif
#              endif // if ( MODEL == HYDRO  &&  !defined SRHD )

//             record the amount of sunk variables removed at the maximum level
               if ( CurrentMaxLv )
               {
                  const real Ekin_old = (real)0.5*( SQR(fluid_old[MOMX]) + SQR(fluid_old[MOMY]) + SQR(fluid_old[MOMZ]) ) / fluid_old[DENS];
                  const real Ekin_new = (real)0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
                  const real Eint_old = fluid_old[ENGY] - Ekin_old - Emag;
                  const real Eint_new = fluid[ENGY]     - Ekin_new - Emag;

#                 pragma omp critical
                  {
                  const int Reset_m1 = Reset-1;
                  CM_BH[Reset_m1].Bondi_SinkMass    += dv *     ( fluid[DENS] - fluid_old[DENS] );
                  CM_BH[Reset_m1].Bondi_SinkMomX    += dv *     ( fluid[MOMX] - fluid_old[MOMX] );
                  CM_BH[Reset_m1].Bondi_SinkMomY    += dv *     ( fluid[MOMY] - fluid_old[MOMY] );
                  CM_BH[Reset_m1].Bondi_SinkMomZ    += dv *     ( fluid[MOMZ] - fluid_old[MOMZ] );
                  CM_BH[Reset_m1].Bondi_SinkMomXAbs += dv * FABS( fluid[MOMX] - fluid_old[MOMX] );
                  CM_BH[Reset_m1].Bondi_SinkMomYAbs += dv * FABS( fluid[MOMY] - fluid_old[MOMY] );
                  CM_BH[Reset_m1].Bondi_SinkMomZAbs += dv * FABS( fluid[MOMZ] - fluid_old[MOMZ] );
                  CM_BH[Reset_m1].Bondi_SinkE       += dv *     ( fluid[ENGY] - fluid_old[ENGY] );
                  CM_BH[Reset_m1].Bondi_SinkEk      += dv *     ( Ekin_new    - Ekin_old        );
                  CM_BH[Reset_m1].Bondi_SinkEt      += dv *     ( Eint_new    - Eint_old        );
                  CM_BH[Reset_m1].Bondi_SinkNCell   ++;
                  }
               } // if ( CurrentMaxLv )

//             store the reset values
               for (int v=0; v<NCOMP_TOTAL; v++)
                  amr->patch[FluSg][lv][PID]->fluid[v][k][j][i] = fluid[v];
            } // if ( Reset != 0 )
         }}} // i,j,k
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // if ( AGN_feedback )

} // FUNCTION : Flu_ResetByUser_API_ClusterMerger



//-------------------------------------------------------------------------------------------------------
// Function    :  GetClusterCenter
// Description :  Get the cluster centers
//
// Note        :  1. Must enable Merger_Coll_LabelCenter
//                2. This routine will only change the cluster centers and velocities if
//                   AdjustPos and AdjustVel (see below) are true and the level is the
//                   maximum level currently refined.
//
// Parameter   :  lv        : Refinement level to search on
//                AdjustPos : If true, adjust the positions of the BHs
//                AdjustVel : If true, adjust the velocities of the BHs
//                BH        : Vector of black holes; BH[c].Pos/Vel are read as the previous center/velocity
//                            and overwritten with the newly computed ones
//
// Return      :  BH[].Pos, BH[].Vel, BH[].Cluster_NPar_close (updated by this function)
//-------------------------------------------------------------------------------------------------------
void GetClusterCenter( const int lv, const bool AdjustPos, const bool AdjustVel, std::vector<BlackHole_t> &BH )
{

// fix the BH position and rest BH
   if ( fixBH )
   {
      for (int d=0; d<3; d++)   BH[0].Pos[d] = amr->BoxCenter[d];
      for (int d=0; d<3; d++)   BH[0].Vel[d] = 0.0;

      return;
   } // if ( fixBH )

   const bool CurrentMaxLv = (  NPatchTotal[lv] > 0  &&  ( lv == MAX_LEVEL || NPatchTotal[lv+1] == 0 )  );
   double pos_min[Merger_Coll_NumBHs][3], DM_Vel[Merger_Coll_NumBHs][3];   // the updated BH position and velocity

// initialize pos_min to be the old center
   for (int c=0; c<Merger_Coll_NumBHs; c++)
      for (int d=0; d<3; d++)    pos_min[c][d] = BH[c].Pos[d];

// Only get the cluster centers and velocities if we are required to adjust them, and if
// this is the maximum level currently refined. Otherwise, this routine skips to the end
// of the function.
   if (  CurrentMaxLv  &&  ( AdjustPos || AdjustVel )  )
   {

      const double dis_exp     = 1e-6*Const_Mpc/UNIT_L;  // to check if the output BH positions of each calculation are close enough
      const int    max_counter = 10;                     // The maximum number of calculations allowed to find the center
      bool converged = false;                            // if the BH positions are close enough, then complete the calculation
      int  counter   = 0;                                // how many times the calculation is performed (minimum: 2, maximum: 10)
      double Cen_new_pre[Merger_Coll_NumBHs][3];
      int    NParClose[Merger_Coll_NumBHs];               // number of particles within 10x the accretion radius of each BH

//###OPTIMIZATION: This loop could be optimized by applying the convergence check to each
//cluster separately to avoid redundant work. Currently, the check is performed on both
//clusters until both are converged.
      while ( converged == false  &&  counter <= max_counter )
      {
         for (int c=0; c<Merger_Coll_NumBHs; c++)
            for (int d=0; d<3; d++)  Cen_new_pre[c][d] = pos_min[c][d];

//       find the particles within "search_factor" times the accretion radius and, on the fly,
//       track each BH's local minimum potential (and the position where it occurs) and the
//       local mass-weighted velocity sums --> no need to gather individual particle data across
//       ranks since amr->Par->Pot[] already holds the true (grid-interpolated) potential
         const double search_factor = 10.0;
         int    num_par[Merger_Coll_NumBHs];         // (each rank) number of particles inside the target region of each cluster
         double Pot_local_min[Merger_Coll_NumBHs];   // (each rank) local minimum potential found so far
         double Pos_at_min[Merger_Coll_NumBHs][3];   // (each rank) position of the particle at Pot_local_min
         double MassSum_local[Merger_Coll_NumBHs];   // (each rank) local sum of particle mass (for the velocity average)
         double VelSum_local[Merger_Coll_NumBHs][3]; // (each rank) local sum of mass*velocity

         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            num_par      [c]    = 0;
            Pot_local_min[c]    = __DBL_MAX__;
            MassSum_local[c]    = 0.0;
            for (int d=0; d<3; d++) { Pos_at_min[c][d] = 0.0;   VelSum_local[c][d] = 0.0; }

            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            {
               const double *EdgeL        = amr->patch[0][lv][PID]->EdgeL;
               const double *EdgeR        = amr->patch[0][lv][PID]->EdgeR;
               const double  patch_pos[3] = { (EdgeL[0]+EdgeR[0])*0.5, (EdgeL[1]+EdgeR[1])*0.5, (EdgeL[2]+EdgeR[2])*0.5 };
               const double  patch_d      = DIST_3D_DBL( EdgeL, EdgeR ) * 0.5;

               if (  DIST_SQR_3D( patch_pos, Cen_new_pre[c] ) > SQR( search_factor*R_acc+patch_d )  )    continue;

               for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               {
                  const long_par ParID         = amr->patch[0][lv][PID]->ParList[p];
                  const real_par ParX_tmp      = amr->Par->PosX[ParID];
                  const real_par ParY_tmp      = amr->Par->PosY[ParID];
                  const real_par ParZ_tmp      = amr->Par->PosZ[ParID];
                  const real_par ParM_tmp      = amr->Par->Mass[ParID];
                  const real_par VelX_tmp      = amr->Par->VelX[ParID];
                  const real_par VelY_tmp      = amr->Par->VelY[ParID];
                  const real_par VelZ_tmp      = amr->Par->VelZ[ParID];
                  const real_par ParPos_tmp[3] = { ParX_tmp, ParY_tmp, ParZ_tmp };

//                Only proceed if the BH belongs to this halo and we are
//                within 10 times the accretion radius
                  if ( CM_ClusterIdx_Cur[amr->Par->AttributeInt[Idx_ParHalo][ParID]] != (long_par)c )    continue;
                  if (  DIST_SQR_3D( ParPos_tmp, Cen_new_pre[c] ) > SQR( search_factor*R_acc )  )        continue;

                  num_par[c] += 1;

                  if ( AdjustPos )
                  {
#                    ifdef STORE_PAR_POT
                     const double Pot_tmp = (double)amr->Par->Pot[ParID];

                     if ( Pot_tmp < Pot_local_min[c] )
                     {
                        Pot_local_min[c] = Pot_tmp;
                        Pos_at_min[c][0] = (double)ParX_tmp;
                        Pos_at_min[c][1] = (double)ParY_tmp;
                        Pos_at_min[c][2] = (double)ParZ_tmp;
                     }
#                    else
                     Aux_Error( ERROR_INFO, "STORE_PAR_POT must be enabled !!\n" );
#                    endif
                  } // if ( AdjustPos )

                  if ( AdjustVel )
                  {
                     MassSum_local[c]    += (double)ParM_tmp;
                     VelSum_local [c][0] += (double)ParM_tmp*(double)VelX_tmp;
                     VelSum_local [c][1] += (double)ParM_tmp*(double)VelY_tmp;
                     VelSum_local [c][2] += (double)ParM_tmp*(double)VelZ_tmp;
                  } // if ( AdjustVel )
               } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
            } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         } // for (int c=0; c<Merger_Coll_NumBHs; c++)

//       collect the number of target particles from each rank
         MPI_Allreduce( num_par, NParClose, Merger_Coll_NumBHs, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

//       find the global potential minimum (and its position) across all ranks
         if ( AdjustPos )
         {
            struct { double val; int rank; } local_minloc[Merger_Coll_NumBHs], global_minloc[Merger_Coll_NumBHs];

            for (int c=0; c<Merger_Coll_NumBHs; c++)
            {
               local_minloc[c].val  = Pot_local_min[c];
               local_minloc[c].rank = MPI_Rank;
            }

            MPI_Allreduce( local_minloc, global_minloc, Merger_Coll_NumBHs, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD );

            for (int c=0; c<Merger_Coll_NumBHs; c++)
            {
               MPI_Bcast( Pos_at_min[c], 3, MPI_DOUBLE, global_minloc[c].rank, MPI_COMM_WORLD );
               for (int d=0; d<3; d++)   pos_min[c][d] = Pos_at_min[c][d];
            }
         } // if ( AdjustPos )

//       calculate the average DM velocity
//       --> use flat 1D buffers (rather than double[Merger_Coll_NumBHs][4]) since the serial-mode
//           MPI_Allreduce substitute (include/Serial.h) just does a flat element-wise copy indexed
//           by Count, which requires SBuf/RBuf to be 1D arrays of that length
         if ( AdjustVel )
         {
            double send[4*Merger_Coll_NumBHs], recv[4*Merger_Coll_NumBHs];

            for (int c=0; c<Merger_Coll_NumBHs; c++)
            {
               send[4*c] = MassSum_local[c];
               for (int d=0; d<3; d++)   send[4*c+1+d] = VelSum_local[c][d];
            }

            MPI_Allreduce( send, recv, 4*Merger_Coll_NumBHs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

            for (int c=0; c<Merger_Coll_NumBHs; c++)
            {
               const double ParM_Tot = recv[4*c];
               for (int d=0; d<3; d++)   DM_Vel[c][d] = recv[4*c+1+d] / ParM_Tot;
            }
         } // if ( AdjustVel )

//       iterate the above calculation until the output BH positions become close enough
         counter += 1;
         double dis[Merger_Coll_NumBHs];

         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            dis[c] = 0.0;
            for (int d=0; d<3; d++)   dis[c] += SQR( pos_min[c][d] - Cen_new_pre[c][d] );
         }

         if ( counter > 1 )
         {
            converged = true;
            for (int c=0; c<Merger_Coll_NumBHs; c++)
               if ( sqrt(dis[c]) >= dis_exp )   converged = false;
         }

         for (int c=0; c<Merger_Coll_NumBHs; c++)   BH[c].Cluster_NPar_close = NParClose[c];
      } // while ( converged == false  &&  counter <= max_counter )
   } // if ( CurrentMaxLv  &&  (AdjustPos  ||  AdjustVel) )

// find the BH particles and adjust their position and velocity
   for (int c=0; c<Merger_Coll_NumBHs; c++)
   {
      double Cen_Tmp[3] = { -__FLT_MAX__, -__FLT_MAX__, -__FLT_MAX__ }; // set to -inf
      double Vel_Tmp[3] = { -__FLT_MAX__, -__FLT_MAX__, -__FLT_MAX__ };

      for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
      {
         if ( amr->Par->Mass[p] < (real_par)0.0 )                                         continue;
         if ( CM_ClusterIdx_Cur[amr->Par->AttributeInt[Idx_ParHalo][p]] != (long_par)c )  continue;
         if ( amr->Par->Type[p] != PTYPE_BLACK_HOLE )                                     continue;

         if ( CurrentMaxLv  &&  AdjustPos )
         {
            amr->Par->PosX[p] = pos_min[c][0];
            amr->Par->PosY[p] = pos_min[c][1];
            amr->Par->PosZ[p] = pos_min[c][2];
         }

         if ( CurrentMaxLv  &&  AdjustVel )
         {
            amr->Par->VelX[p] = DM_Vel[c][0];
            amr->Par->VelY[p] = DM_Vel[c][1];
            amr->Par->VelZ[p] = DM_Vel[c][2];
         }

         Cen_Tmp[0] = amr->Par->PosX[p];
         Cen_Tmp[1] = amr->Par->PosY[p];
         Cen_Tmp[2] = amr->Par->PosZ[p];
         Vel_Tmp[0] = amr->Par->VelX[p];
         Vel_Tmp[1] = amr->Par->VelY[p];
         Vel_Tmp[2] = amr->Par->VelZ[p];

         break;
      } // for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)

//    use MPI_MAX since Cen_Tmp[] and Vel_Tmp[] are initialized as -inf
      MPI_Allreduce( Cen_Tmp, BH[c].Pos, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
      MPI_Allreduce( Vel_Tmp, BH[c].Vel, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
   } // for (int c=0; c<Merger_Coll_NumBHs; c++)

   if ( CurrentMaxLv  &&  AdjustPos )
   {
      const bool TimingSendPar_No = false;
      Par_PassParticle2Sibling( lv, TimingSendPar_No );
      Par_PassParticle2Son_MultiPatch( lv, PAR_PASS2SON_EVOLVE, TimingSendPar_No, NULL_INT, NULL );
   }

} // FUNCTION : GetClusterCenter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetJetDirection
// Description :  Set the jet direction
//
// Parameter   :  TimeNew : Current physical time (system has been updated from TimeOld to TimeNew in EvolveLevel())
//
// Return      :  CM_BH[].Jet_Vec
//-------------------------------------------------------------------------------------------------------
void SetJetDirection( const double TimeNew, const int lv, const int FluSg )
{

   switch ( JetDirection_case )
   {
      case 1: // fixed at x-axis
         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            CM_BH[c].Jet_Vec[0] = 1.0;
            CM_BH[c].Jet_Vec[1] = 0.0;
            CM_BH[c].Jet_Vec[2] = 0.0;
         }
         break;

      case 2: // import from table
         {
             const double Time_period      = CM_Jet_Time_table[JetDirection_NBin-1];
             const double Time_interpolate = fmod( TimeNew, Time_period );
             for (int c=0; c<Merger_Coll_NumBHs; c++)
             {
                const double theta = Mis_InterpolateFromTable( JetDirection_NBin, CM_Jet_Time_table, CM_Jet_Theta_table[c], Time_interpolate );
                const double phi   = Mis_InterpolateFromTable( JetDirection_NBin, CM_Jet_Time_table, CM_Jet_Phi_table[c],   Time_interpolate );

                CM_BH[c].Jet_Vec[0] = cos(theta);
                CM_BH[c].Jet_Vec[1] = sin(theta)*cos(phi);
                CM_BH[c].Jet_Vec[2] = sin(theta)*sin(phi);
             }
         }
         break;

      case 3: // align with angular momentum
         {
//          angular momentum inside the accretion radius, per rank and total
            const double dh = amr->dh[lv];
            const double dv = CUBE(dh);
            double ang_mom    [Merger_Coll_NumBHs][3];
            double ang_mom_sum[Merger_Coll_NumBHs][3];

            for (int c=0; c<Merger_Coll_NumBHs; c++)
               for (int d=0; d<3; d++)
                  ang_mom[c][d] = 0.0;

            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            {
               for (int k=0; k<PS1; k++)
               for (int j=0; j<PS1; j++)
               for (int i=0; i<PS1; i++)
               {
                  const double pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (0.5+i)*dh,
                                          amr->patch[0][lv][PID]->EdgeL[1] + (0.5+j)*dh,
                                          amr->patch[0][lv][PID]->EdgeL[2] + (0.5+k)*dh };

                  for (int c=0; c<Merger_Coll_NumBHs; c++)
                  {
                     if (  DIST_SQR_3D( pos, CM_BH[c].Pos ) <= SQR( R_acc )  )
                     {
                        double dr[3] = { pos[0]-CM_BH[c].Pos[0], pos[1]-CM_BH[c].Pos[1], pos[2]-CM_BH[c].Pos[2] };
                        ang_mom[c][0] += dv * ( dr[1]*amr->patch[FluSg][lv][PID]->fluid[MOMZ][k][j][i] - dr[2]*amr->patch[FluSg][lv][PID]->fluid[MOMY][k][j][i] );
                        ang_mom[c][1] += dv * ( dr[2]*amr->patch[FluSg][lv][PID]->fluid[MOMX][k][j][i] - dr[0]*amr->patch[FluSg][lv][PID]->fluid[MOMZ][k][j][i] );
                        ang_mom[c][2] += dv * ( dr[0]*amr->patch[FluSg][lv][PID]->fluid[MOMY][k][j][i] - dr[1]*amr->patch[FluSg][lv][PID]->fluid[MOMX][k][j][i] );
                     }
                  } // for (int c=0; c<Merger_Coll_NumBHs; c++)
               } // for (int k=0; k<PS1; k++); for (int j=0; j<PS1; j++); for (int i=0; i<PS1; i++)
            } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

            for (int c=0; c<Merger_Coll_NumBHs; c++)
            {
               MPI_Allreduce( ang_mom[c], ang_mom_sum[c], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

               const double ang_mom_norm = sqrt( SQR(ang_mom_sum[c][0]) + SQR(ang_mom_sum[c][1]) + SQR(ang_mom_sum[c][2]) );

               for (int d=0; d<3; d++)   CM_BH[c].Jet_Vec[d] = ang_mom_sum[c][d] / ang_mom_norm;
            } // for (int c=0; c<Merger_Coll_NumBHs; c++)
         }
         break;

      default:
         Aux_Error( ERROR_INFO, "Unsupported JetDirection_case %d [1/2/3] !!\n", JetDirection_case );
         break;
   } // switch ( JetDirection_case )

} // FUNCTION : SetJetDirection



#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY  &&  defined MASSIVE_PARTICLES )
