#include "GAMER.h"
extern FieldIdx_t ParTypeIdx_DiskHeating;
#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFile
// Description :  Initialize particle attributes from a file
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box even for the periodic BC
//                3. Particles lying outside the active region will be removed later by Par_Aux_InitCheck()
//                   if non-periodic B.C. is adopted
//                4. Particles loaded here are only temporarily stored in this rank
//                   --> They will be redistributed when calling Par_FindHomePatch_UniformGrid()
//                       and LB_Init_LoadBalance()
//                   --> So there is no constraint on which particles should be set by this function
//                5. Currently the target file name is fixed to "PAR_IC"
//                6. The data format of the PAR_IC file is controlled by the runtime parameter "PAR_IC_FORMAT"
//                   --> PAR_IC_FORMAT_ATT_ID: [particle attribute][particle id] in a row-major order
//                       PAR_IC_FORMAT_ID_ATT: [particle id][particle attribute] in a row-major order
//                7  Currently it only loads particle mass, position x/y/z, and velocity x/y/z
//                   (and must be in the same order of PAR_MASS, PAR_POSX/Y/Z, and PAR_VELX/Y/Z)
//                   --> The mass of all particles can be set to PAR_IC_MASS instead (by having PAR_IC_MASS>=0.0)
//                       --> In this case, the PAR_IC file should exclude the partice mass data
//                8. For LOAD_BALANCE, the number of particles in each rank must be set in advance
//                   --> Currently it's set by Init_Parallelization()
//
// Parameter   :  None
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, AllAttribute
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_DiskHeating( const long NPar_ThisRank, const long NPar_AllRank,
                                real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                real *ParType, real *AllAttribute[PAR_NATT_TOTAL])
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check
//   if ( amr->Par->ParICFormat != PAR_IC_FORMAT_ID_ATT  &&  amr->Par->ParICFormat != PAR_IC_FORMAT_ATT_ID )
//      Aux_Error( ERROR_INFO, "unknown data format in PAR_IC (%d) !!\n", amr->Par->ParICFormat );

   const char FileName[]     = "PAR_IC";
   const long NParAllRank    = amr->Par->NPar_Active_AllRank;
         long NParThisRank   = amr->Par->NPar_AcPlusInac;            // cannot be "const" due to MPI_Allgather()
   const int  NParAtt        = 8;                                    // mass, pos*3, vel*3, type
   const int  NParAttPerLoad = 1;


// check
   if ( !Aux_CheckFileExist(FileName) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

   FILE *FileTemp = fopen( FileName, "rb" );

   fseek( FileTemp, 0, SEEK_END );

   const long ExpectSize = long(NParAtt)*NParAllRank*sizeof(real);
   const long FileSize   = ftell( FileTemp );
   if ( FileSize != ExpectSize )
      Aux_Error( ERROR_INFO, "size of the file <%s> = %ld != expect = %ld !!\n",
                 FileName, FileSize, ExpectSize );

   fclose( FileTemp );

   MPI_Barrier( MPI_COMM_WORLD );


// set the file offset for this rank
   long NPar_EachRank[MPI_NRank], NPar_Check=0, FileOffset=0;

   MPI_Allgather( &NParThisRank, 1, MPI_LONG, NPar_EachRank, 1, MPI_LONG, MPI_COMM_WORLD );

// check if the total number of particles is correct
   for (int r=0; r<MPI_NRank; r++)  NPar_Check += NPar_EachRank[r];
   if ( NPar_Check != NParAllRank )
      Aux_Error( ERROR_INFO, "total number of particles found (%ld) != expect (%ld) !!\n", NPar_Check, NParAllRank );

   for (int r=0; r<MPI_Rank; r++)   FileOffset = FileOffset + long(NParAttPerLoad)*NPar_EachRank[r]*sizeof(real);


// load data
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading data ... " );

   real *ParData_ThisRank = new real [ NParThisRank*NParAtt ];

// note that fread() may fail for large files if sizeof(size_t) == 4 instead of 8
   FILE *File = fopen( FileName, "rb" );

   for (int v=0; v<NParAtt; v+=NParAttPerLoad)
   {
      fseek( File, FileOffset+v*NParAllRank*sizeof(real), SEEK_SET );
      fread( ParData_ThisRank+v*NParThisRank, sizeof(real), long(NParAttPerLoad)*NParThisRank, File );
   }

   fclose( File );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// store data into the particle repository
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Storing data into particle repository ... " );

   real *ParData1 = new real [NParAtt];

   for (long p=0; p<NParThisRank; p++)
   {
//    collect data for the target particle
//    [att][id]
      for (int v=0; v<NParAtt; v++)
         ParData1[v] = ParData_ThisRank[ v*NParThisRank + p ];


//    assuming that the orders of the particle attributes stored on the disk and in Par->Attribute[] are the same

      ParMass[p] = ParData1[0];
      ParPosX[p] = ParData1[1];
      ParPosY[p] = ParData1[2];
      ParPosZ[p] = ParData1[3];
      ParVelX[p] = ParData1[4];
      ParVelY[p] = ParData1[5];
      ParVelZ[p] = ParData1[6];
      ParType[p] = ParData1[7];
      AllAttribute[ParTypeIdx_DiskHeating][p] = ParData1[7];

//    synchronize all particles to the physical time at the base level
      amr->Par->Time[p] = Time[0];
   }

   delete [] ParData_ThisRank;
   delete [] ParData1;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_DiskHeating



#endif // #ifdef PARTICLE
