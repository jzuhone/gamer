#include "GAMER.h"
#include "BlackHole.h"




// problem-specific global variables
// =======================================================================================
extern int      Merger_Coll_NumBHs;
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_ClusterMerger
// Description :  Record the cluster centers
//
// Note        :  1. Invoked by main() using the function pointer "Aux_Record_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__RECORD_USER"
//                3. This function will be called both during the program initialization and after each full update
//                4. Must enable Merger_Coll_LabelCenter
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_ClusterMerger()
{

   const char FileName[] = "Record__ClusterCenter";
   static bool FirstTime = true;

// header
   if ( FirstTime )
   {
      if ( MPI_Rank == 0 )
      {
         if (  Aux_CheckFileExist( FileName )  )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         FILE *File_User = fopen( FileName, "a" );
         fprintf( File_User, "# ClusterCen_x/y/z : BH position [code length]\n" );
         fprintf( File_User, "# BHVel_?          : BH velocity [km/s]\n" );
         fprintf( File_User, "# GasVel_?         : Gas velocity [km/s]\n" );
         fprintf( File_User, "# RelativeVel      : Relative velocity between BH and gas for each cluster inside the accretion radius [km/s]\n" );
         fprintf( File_User, "# SoundSpeed       : Average sound speed inside the accretion radius [km/s]\n" );
         fprintf( File_User, "# GasDens          : Average gas density inside the accretion radius [cgs]\n" );
         fprintf( File_User, "# mass_BH          : BH mass [MSun]\n" );
         fprintf( File_User, "# Mdot_?_BH        : BH accretion rate [cgs]\n" );
         fprintf( File_User, "# NVoidCell        : Total number of finest cells within the feedback region\n" );
         fprintf( File_User, "# Mom?Inj          : Total momentum change in the feedback region [cgs]\n" );
         fprintf( File_User, "# Mom?InjAbs       : Total |momentum| change in the feedback region [cgs]\n" );
         fprintf( File_User, "# EInj_exp         : The expected amount of injected energy in the feedback region [erg]\n" );
         fprintf( File_User, "# E_Inj            : Total injected energy in the feedback region [erg]\n" );
         fprintf( File_User, "# E_Inj_err        : The relative error of the total injected energy\n" );
         fprintf( File_User, "# Ek_Inj           : Total injected kinetic energy in the feedback region [erg]\n" );
         fprintf( File_User, "# Et_Inj           : Total injected thermal energy in the feedback region [erg]\n" );
         fprintf( File_User, "# PowerInj         : Total injected power in the feedback region [cgs]\n" );
         fprintf( File_User, "# MInjexp          : The expected amount of injected gas mass in the feedback region [Msun]\n" );
         fprintf( File_User, "# MassInj          : Total mass change in the feedback region [Msun]\n" );
         fprintf( File_User, "# M_Inj_err        : The relative error of the total injected mass\n" );
         fprintf( File_User, "# Mdot_Inj         : The feedback injection rate of mass [cgs]\n" );
         fprintf( File_User, "# Pdot_Inj         : The feedback injection rate of momentum [cgs]\n" );
         fprintf( File_User, "# Edot_Inj         : The feedback injection rate of total energy [cgs]\n" );
         fprintf( File_User, "# Jet_Vec_?        : Jet inject direction\n" );
         fprintf( File_User, "# num_par_sum      : Total number of particles inside the accretion region\n" );
         fprintf( File_User, "# ColdGasMass      : Cold gas mass inside the accretion radius [Msun]\n" );
         fprintf( File_User, "#\n" );

         fprintf( File_User, "#%23s %24s", "[  1]", "[  2]");
         const int N_att = 40;
         int idx = 3;
         for (int c=0; c<Merger_Coll_NumBHs; c++)
         for (int i=0; i<N_att; i++)
         {
            fprintf( File_User, " %19s[%3d]", "", idx );
            idx += 1;
         }
         fprintf( File_User, "\n" );
         fprintf( File_User, "#%23s %24s",  "Time", "Step" );
         for (int c=0; c<Merger_Coll_NumBHs; c++)
         {
            fprintf( File_User, " %23s%1d %23s%1d %23s%1d", "ClusterCen_x", c, "ClusterCen_y", c, "ClusterCen_z", c );
            fprintf( File_User, " %23s%1d %23s%1d %23s%1d", "BHVel_x",      c, "BHVel_y",      c, "BHVel_z",      c );
            fprintf( File_User, " %23s%1d %23s%1d %23s%1d", "GasVel_x",     c, "GasVel_y",     c, "GasVel_z",     c );
            fprintf( File_User, " %23s%1d %23s%1d %23s%1d", "RelativeVel",  c, "SoundSpeed",   c, "GasDens",      c );
            fprintf( File_User, " %23s%1d",                 "mass_BH",      c );
            fprintf( File_User, " %23s%1d %23s%1d %23s%1d", "Mdot_tot_BH",  c, "Mdot_hot_BH",  c, "Mdot_cold_BH", c );
            fprintf( File_User, " %23s%1d",                 "NVoidCell",    c );
            fprintf( File_User, " %23s%1d %23s%1d %23s%1d", "MomXInj",      c, "MomYInj",      c, "MomZInj",      c );
            fprintf( File_User, " %23s%1d %23s%1d %23s%1d", "MomXInjAbs",   c, "MomYInjAbs",   c, "MomZInjAbs",   c );
            fprintf( File_User, " %23s%1d %23s%1d %23s%1d", "EInj_exp",     c, "E_Inj",        c, "E_Inj_err",    c );
            fprintf( File_User, " %23s%1d %23s%1d %23s%1d", "Ek_Inj",       c, "Et_Inj",       c, "PowerInj",     c );
            fprintf( File_User, " %23s%1d %23s%1d %23s%1d", "MInjexp",      c, "MassInj",      c, "M_Inj_err",    c );
            fprintf( File_User, " %23s%1d %23s%1d %23s%1d", "Mdot_Inj",     c, "Pdot_Inj",     c, "Edot_Inj",     c );
            fprintf( File_User, " %23s%1d %23s%1d %23s%1d", "Jet_Vec_x",    c, "Jet_Vec_y",    c, "Jet_Vec_z",    c );
            fprintf( File_User, " %23s%1d %23s%1d",         "num_par_sum",  c, "ColdGasMass",  c );
         }
         fprintf( File_User, "\n" );
         fclose( File_User );
      } // if ( MPI_Rank == 0 )

      FirstTime = false;
   } // if ( FirstTime )


// sum over the variables and convert units
   int    SinkNCell_Sum[Merger_Coll_NumBHs];
   double Mass_Sum[Merger_Coll_NumBHs];
   double MomX_Sum[Merger_Coll_NumBHs], MomY_Sum[Merger_Coll_NumBHs], MomZ_Sum[Merger_Coll_NumBHs];
   double MomXAbs_Sum[Merger_Coll_NumBHs], MomYAbs_Sum[Merger_Coll_NumBHs], MomZAbs_Sum[Merger_Coll_NumBHs];
   double E_Sum[Merger_Coll_NumBHs], Ek_Sum[Merger_Coll_NumBHs], Et_Sum[Merger_Coll_NumBHs];
   double E_power_inj[Merger_Coll_NumBHs];

   for (int c=0; c<Merger_Coll_NumBHs; c++)
   {
      MPI_Reduce( &CM_BH[c].Bondi_SinkNCell,   &SinkNCell_Sum[c], 1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_BH[c].Bondi_SinkMass,    &Mass_Sum[c],      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_BH[c].Bondi_SinkMomX,    &MomX_Sum[c],      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_BH[c].Bondi_SinkMomY,    &MomY_Sum[c],      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_BH[c].Bondi_SinkMomZ,    &MomZ_Sum[c],      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_BH[c].Bondi_SinkMomXAbs, &MomXAbs_Sum[c],   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_BH[c].Bondi_SinkMomYAbs, &MomYAbs_Sum[c],   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_BH[c].Bondi_SinkMomZAbs, &MomZAbs_Sum[c],   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_BH[c].Bondi_SinkE,       &E_Sum[c],         1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_BH[c].Bondi_SinkEk,      &Ek_Sum[c],        1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_BH[c].Bondi_SinkEt,      &Et_Sum[c],        1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

      Mass_Sum   [c] *= UNIT_M / Const_Msun;
      MomX_Sum   [c] *= UNIT_M * UNIT_V;
      MomY_Sum   [c] *= UNIT_M * UNIT_V;
      MomZ_Sum   [c] *= UNIT_M * UNIT_V;
      MomXAbs_Sum[c] *= UNIT_M * UNIT_V;
      MomYAbs_Sum[c] *= UNIT_M * UNIT_V;
      MomZAbs_Sum[c] *= UNIT_M * UNIT_V;
      E_Sum      [c] *= UNIT_E;
      Ek_Sum     [c] *= UNIT_E;
      Et_Sum     [c] *= UNIT_E;
      CM_BH[c].E_inj_exp *= UNIT_E;
      CM_BH[c].M_inj_exp *= UNIT_M / Const_Msun;
   } // for (int c=0; c<Merger_Coll_NumBHs; c++)

   for (int c=0; c<Merger_Coll_NumBHs; c++)   E_power_inj[c] = E_Sum[c]/(dTime_Base*UNIT_T);


// output the properties of the cluster centers
   if ( MPI_Rank == 0 )
   {
      FILE *File_User = fopen( FileName, "a" );
      fprintf( File_User, "%24.16e %24ld", Time[0], Step );
      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         for (int d=0; d<3; d++)   fprintf( File_User, " %24.16e", CM_BH[c].Pos[d] );
         for (int d=0; d<3; d++)   fprintf( File_User, " %24.16e", CM_BH[c].Vel[d]*UNIT_V/(Const_km/Const_s) );
         for (int d=0; d<3; d++)   fprintf( File_User, " %24.16e", CM_BH[c].RAcc_GasVel[d]*UNIT_V/(Const_km/Const_s) );
         fprintf( File_User, " %24.16e %24.14e %24.16e", CM_BH[c].RAcc_RelativeVel*UNIT_V/(Const_km/Const_s), CM_BH[c].RAcc_SoundSpeed*UNIT_V/(Const_km/Const_s), CM_BH[c].RAcc_GasDens*UNIT_D );
         fprintf( File_User, " %24.16e",                 CM_BH[c].Mass*UNIT_M/Const_Msun );
         fprintf( File_User, " %24.16e %24.14e %24.16e", CM_BH[c].Mdot_tot*UNIT_M/UNIT_T, CM_BH[c].Mdot_hot*UNIT_M/UNIT_T, CM_BH[c].Mdot_cold*UNIT_M/UNIT_T );
         fprintf( File_User, " %24d",                    SinkNCell_Sum[c] );
         fprintf( File_User, " %24.16e %24.16e %24.16e", MomX_Sum[c], MomY_Sum[c], MomZ_Sum[c] );
         fprintf( File_User, " %24.16e %24.16e %24.16e", MomXAbs_Sum[c], MomYAbs_Sum[c], MomZAbs_Sum[c] );
         fprintf( File_User, " %24.16e %24.16e %24.16e", CM_BH[c].E_inj_exp, E_Sum[c], (E_Sum[c]-CM_BH[c].E_inj_exp)/CM_BH[c].E_inj_exp );
         fprintf( File_User, " %24.16e %24.16e %24.16e", Ek_Sum[c], Et_Sum[c], E_power_inj[c] );
         fprintf( File_User, " %24.16e %24.16e %24.16e", CM_BH[c].M_inj_exp, Mass_Sum[c], (Mass_Sum[c]-CM_BH[c].M_inj_exp)/CM_BH[c].M_inj_exp );
         fprintf( File_User, " %24.16e %24.16e %24.16e", CM_BH[c].Jet_Mdot*UNIT_M/UNIT_T, CM_BH[c].Jet_Pdot*UNIT_M*UNIT_V/UNIT_T,
                                                         CM_BH[c].Jet_Edot*UNIT_E/UNIT_T );
         for (int d=0; d<3; d++)   fprintf( File_User, " %24.16e", CM_BH[c].Jet_Vec[d] );
         fprintf( File_User, " %24d %24.16e",            CM_BH[c].Cluster_NPar_close, CM_BH[c].RAcc_ColdGasMass*UNIT_M/Const_Msun );
      }
      fprintf( File_User, "\n" );
      fclose( File_User );
   } // if ( MPI_Rank == 0 )


// reset the cumulative variables to zero
   for (int c=0; c<Merger_Coll_NumBHs; c++)   CM_BH[c].ResetRecordDiagnostics();

} // FUNCTION : Aux_Record_ClusterMerger
