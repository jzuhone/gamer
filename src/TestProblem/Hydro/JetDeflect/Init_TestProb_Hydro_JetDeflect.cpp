#include "GAMER.h"



static FieldIdx_t JetFieldIdx  = Idx_Undefined;
static FieldIdx_t ICMFieldIdx  = Idx_Undefined;
static FieldIdx_t LobeFieldIdx = Idx_Undefined;
static FieldIdx_t IntFieldIdx  = Idx_Undefined;

// problem-specific global variables
// =======================================================================================
// background parameters
static double   ICM_Density;             // ICM density
static double   Jump_Position_x;         // position of interface
static double   Jump_Position_y;         // position of interface
static double   Jump_Angle;              // inclination of interface
static double   Jump_Tangent;            // tangent of interface angle
static double   Jump_Width;              // width of jump
static double   Amb_Pressure;            // ambient pressure
static double   Lobe_ICM_Ratio;          // ratio of lobe/ICM densities
static double   Lobe_Density;            // lobe density
static bool     Wall_Off;                // turn the ICM wall off
static bool     ICM_Blob;                // turn the ICM blob on
static double   Blob_PosX;               // x-position of blob
static double	Blob_PosY;    	      	 // y-position of blob
static double   Blob_PosZ;               // z-position of blob
static double   Blob_Radius;             // radius of blob
static double   Blob_VelY;               // y-velocity of blob
static double   Jump_Sine;               // sine of jump angle
static double   Jump_Cosine;             // cosine of jump angle

// jet parameters
static bool     Jet_Fire;                // [true/false]: jet on/off
static int      Jet_Axis;                // axis of jet
static double   Jet_Velocity;            // jet velocity (units of c)
static double   Jet_VelSlope;            // Slope of velocity gradient across jet
static double   Jet_VelCenter;           // jet central velocity
static double   Jet_Radius;              // radius of jet
static double   Jet_Position;            // position of jet
static double   Jet_Lobe_Ratio;          // ratio of jet/lobe densities
static double   Jet_Center[2];           // jet central coordinates
static double   Jet_Density;             // jet density
static double   Jet_Gamma;               // jet relativistic gamma
static double   Jet_PrecessAngle;        // jet precession angle
static double   Jet_PrecessPeriod;       // jet precession period
static double   Jet_PrecessOmega;        // jet precession omega
static double   Jet_Cosine;              // jet cosine
static double   Jet_Sine;                // jet sine
static double   Jump_Sine;
static double   Jump_Cosine;

// =======================================================================================

static void JetBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                   const int GhostSize, const int idx[], const double pos[], const double Time,
                   const int lv, const int TFluVarIdxList[], double AuxArray[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


// errors
#  ifndef SRHD
   Aux_Error( ERROR_INFO, "SRHD must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );

   if ( ICM_Blob && !Wall_Off )
      Aux_Error( ERROR_INFO, "To use ICM_Blob, Wall_Off must be set !!\n");
		
// warnings
   if ( MPI_Rank == 0 )
   {
      for (int s=0; s<6; s++)
         if ( OPT__BC_FLU[s] != BC_FLU_OUTFLOW )
            Aux_Message( stderr, "WARNING : it's recommended to use the outflow BC (currently OPT__BC_FLU[%d] = %d != 2)\n",
                         s, OPT__BC_FLU[s] );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",    &VARIABLE,          DEFAULT,       MIN,            MAX          );
// ************************************************************************************************************************

// background parameters
   ReadPara->Add( "ICM_Density",        &ICM_Density,       NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jump_Position_x",    &Jump_Position_x,   NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jump_Position_y",    &Jump_Position_y,   NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jump_Angle",         &Jump_Angle,        NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Amb_Pressure",       &Amb_Pressure,      NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Lobe_ICM_Ratio",     &Lobe_ICM_Ratio,    NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jump_Width",         &Jump_Width,        NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Wall_Off",           &Wall_Off,          false,         Useless_bool,   Useless_bool );
   ReadPara->Add( "ICM_Blob",           &ICM_Blob,          false,         Useless_bool,   Useless_bool );
   ReadPara->Add( "Blob_PosX",          &Blob_PosX,         NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Blob_PosY",          &Blob_PosY,         NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Blob_Radius",        &Blob_Radius,       NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Blob_VelY",          &Blob_VelY,         NoDef_double,  NoMin_double,   NoMax_double );
   
// jet parameters
   ReadPara->Add( "Jet_Fire",           &Jet_Fire,          false,         Useless_bool,   Useless_bool );
   ReadPara->Add( "Jet_Axis",           &Jet_Axis,          1,             0,              1            );
   ReadPara->Add( "Jet_Lobe_Ratio",     &Jet_Lobe_Ratio,    NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jet_Radius",         &Jet_Radius,        NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jet_Position",       &Jet_Position,      NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jet_Velocity",       &Jet_Velocity,      NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jet_VelSlope",       &Jet_VelSlope,      0.0,           NoMin_double,   0.0          );
   ReadPara->Add( "Jet_VelCenter",      &Jet_VelCenter,     NoDef_double,  0.0,            NoMax_double );
   ReadPara->Add( "Jet_PrecessAngle",   &Jet_PrecessAngle,  0.0,           0.0,            NoMax_double );
   ReadPara->Add( "Jet_PrecessPeriod",  &Jet_PrecessPeriod, 0.0,           0.0,            NoMax_double );

   ReadPara->Read( FileName );

   delete ReadPara;

   Jet_Gamma = sqrt(1.0+SQR(Jet_Velocity));

// (1-2) convert to code unit
   Jet_Velocity      *= Const_c   / UNIT_V;
   ICM_Density       *= 1.0       / UNIT_D;
   Jet_Radius        *= Const_kpc / UNIT_L;
   Blob_PosX         *= Const_kpc / UNIT_L;
   Blob_PosY  	     *=	Const_kpc / UNIT_L;
   Blob_Radius       *= Const_kpc / UNIT_L;
   Blob_VelY         *= Const_km  / UNIT_V;
   Jet_Position      *= Const_kpc / UNIT_L;
   Jump_Position_x   *= Const_kpc / UNIT_L;
   Jump_Position_y   *= Const_kpc / UNIT_L;
   Amb_Pressure      *= 1.0       / UNIT_P;
   Jump_Width        *= Const_kpc / UNIT_L;
   Jet_PrecessPeriod *= 1000.0*Const_yr / UNIT_T;


// (2) set the problem-specific derived parameters
   Jump_Tangent     = tan( Jump_Angle*M_PI/180.0 );
   Jump_Sine        = sin( Jump_Angle*M_PI/180.0 );
   Jump_Cosine      = cos( Jump_Angle*M_PI/180.0 );
   Lobe_Density     = Lobe_ICM_Ratio*ICM_Density;
   Jet_Density      = Jet_Lobe_Ratio*Lobe_Density;
   Jet_Center[0]    = Jet_Position;
   Jet_Center[1]    = amr->BoxCenter[2];
   Jet_PrecessOmega = Jet_PrecessPeriod > 0.0 ? 2.0*M_PI/Jet_PrecessPeriod : 0.0;
   Jet_Cosine       = cos( Jet_PrecessAngle*M_PI/180.0 );
   Jet_Sine         = sin( Jet_PrecessAngle*M_PI/180.0 );
   Blob_PosZ        = amr->BoxCenter[2];
  

// (3) reset other general-purpose parameters
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 100.0*Const_kyr / UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
     Aux_Message( stdout, "=============================================================================\n " );
     Aux_Message( stdout, "  test problem ID = %d\n",              TESTPROB_ID                      );
     Aux_Message( stdout, "  ICM_Density     = %14.7e g/cm^3\n",   ICM_Density*UNIT_D               );
     Aux_Message( stdout, "  Jump_Position_x = %14.7e kpc\n",      Jump_Position_x*UNIT_L/Const_kpc );
     Aux_Message( stdout, "  Jump_Position_y = %14.7e kpc\n",      Jump_Position_y*UNIT_L/Const_kpc );
     Aux_Message( stdout, "  Jump_Angle      = %14.7e degree\n",   Jump_Angle                       );
     Aux_Message( stdout, "  Jump_Width      = %14.7e kpc\n",      Jump_Width*UNIT_L/Const_kpc      );
     Aux_Message( stdout, "  Amb_Pressure    = %14.7e erg/cm^3\n", Amb_Pressure*UNIT_P              );
     Aux_Message( stdout, "  Lobe_ICM_Ratio  = %14.7e\n",          Lobe_ICM_Ratio                   );
     Aux_Message( stdout, "  Lobe_Density    = %14.7e g/cm^3\n",   Lobe_Density*UNIT_D              );
   }

   if ( Jet_Fire && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_Radius      = %14.7e kpc\n",      Jet_Radius*UNIT_L/Const_kpc      );
     Aux_Message( stdout, "  Jet_Position    = %14.7e kpc\n",      Jet_Position*UNIT_L/Const_kpc    );
     Aux_Message( stdout, "  Jet_Velocity    = %14.7e c\n",        Jet_Velocity*UNIT_V/Const_c      );
     Aux_Message( stdout, "  Jet_VelSlope    = %14.7e\n",          Jet_VelSlope                     );
     Aux_Message( stdout, "  Jet_Density     = %14.7e g/cm^3\n",   Jet_Density*UNIT_D               );
     Aux_Message( stdout, "  Jet_Lobe_Ratio  = %14.7e\n",          Jet_Lobe_Ratio                   );
   }

   if ( MPI_Rank == 0 )
   {
     Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   double d, uy = 0.0;
  
   real PriReal[NCOMP_TOTAL];
   
   if ( Wall_Off ) {

      d = Lobe_Density;
         
      if ( ICM_Blob ) {

	 double rr = sqrt( SQR(x-Blob_PosX) + SQR(y-Blob_PosY) + SQR(z-Blob_PosZ) );

	 if ( rr <= Blob_Radius ) {

	    d = ICM_Density;
	    uy = Blob_VelY;
	    
	 }
      }
      
   } else {
     
      double xx = (x - Jump_Position_x)*Jump_Cosine - (y - Jump_Position_y)*Jump_Sine;
      double xw = xx/Jump_Width;
      
      if (xw > 200.0)
         d = Lobe_Density;
      else
         d = (ICM_Density + Lobe_Density*exp(xw)) / (1.0 + exp(xw));

   }
   
   PriReal[0] = (real)d;
   PriReal[1] = 0.0;
   PriReal[2] = (real)uy;
   PriReal[3] = 0.0;
   PriReal[4] = (real)Amb_Pressure;

   Hydro_Pri2Con( PriReal, fluid, false, PassiveNorm_NVar, PassiveNorm_VarIdx,
                  EoS_DensPres2Eint_CPUPtr, EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

   double ICM_x, lobe_x;
   
   if ( Wall_Off ) {

     ICM_x = 0.0;
     lobe_x = 1.0;
     
   } else {
     
     ICM_x = (d-0.4*ICM_Density)/(0.4*ICM_Density);
     lobe_x = (2.5*Lobe_Density-d)/(1.25*Lobe_Density);
     if ( ICM_x  < 0.0 ) ICM_x  = 0.0;
     if ( lobe_x < 0.0 ) lobe_x = 0.0;
     if ( ICM_x  > 1.0 ) ICM_x  = 1.0;
     if ( lobe_x > 1.0 ) lobe_x = 1.0;

   }
   
   fluid[JetFieldIdx ] = 0.0;
   fluid[ICMFieldIdx ] = (real)(ICM_x*d);
   fluid[LobeFieldIdx] = (real)(lobe_x*d);
   fluid[IntFieldIdx ] = (real)((1-ICM_x-lobe_x)*d);

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  JetBC
// Description :  Set the external boundary condition for the JetDeflect problem. On the -y
//                boundary, inside the jet region, the jet inflows with the jet density and a
//                velocity profile that peaks in the jet center and linearly falls off to the
//                edge of the jet. The jet passive scalar density is set to the jet density
//                within this region. Outside the jet region, the boundary condition is "diode",
//                meaning outflow-only.
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//
// Parameter   :  Array          : Array to store the prepared data including ghost zones
//                ArraySize      : Size of Array including the ghost zones on each side
//                fluid          : Fluid fields to be set
//                NVar_Flu       : Number of fluid variables to be prepared
//                GhostSize      : Number of ghost zones
//                idx            : Array indices
//                pos            : Physical coordinates
//                Time           : Physical time
//                lv             : Refinement level
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                AuxArray       : Auxiliary array
//
// Return                        : fluid
//-------------------------------------------------------------------------------------------------------
void JetBC( real Array[], const int ArraySize[], real BVal[], const int NVar_Flu,
            const int GhostSize, const int idx[], const double pos[], const double Time,
            const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

   int i, j, k, i_ref, j_ref, mom_axis, mom_x, mom_y, jet_x, jet_y;
   
   i = idx[0];
   j = idx[1];
   k = idx[2];
      
   switch ( Jet_Axis ) {
      case 0:
	 mom_axis = MOMX;
	 mom_x    = MOMY;
	 mom_y    = MOMZ;
	 jet_x = 1;
	 jet_y = 2;
	 i_ref = GhostSize;
	 j_ref = j;
         break;
      case 1:
	 mom_axis = MOMY;
         mom_x    = MOMX;
         mom_y    = MOMZ;
	 jet_x = 0;
	 jet_y = 2;
	 i_ref = i;
	 j_ref = GhostSize;
	 break;
   }
   
   real PriReal[NCOMP_TOTAL];

   int TFluVarIdx;

   double rad = sqrt( SQR(pos[jet_x]-Jet_Center[0]) + SQR(pos[jet_y]-Jet_Center[1]) );

// 1D array -> 3D array
   typedef real (*vla)[ ArraySize[2] ][ ArraySize[1] ][ ArraySize[0] ];
   vla Array3D = ( vla )Array;

   double x = rad/Jet_Radius;

   if ( Jet_Fire  &&  x <= 1.0 )
   {
      const double u_jet      = Jet_Velocity*( Jet_VelSlope*x+Jet_VelCenter );
      const double cos_phi    = cos( Jet_PrecessOmega*Time );
      const double sin_phi    = sin( Jet_PrecessOmega*Time );
      const double LntzFact   = sqrt( 1.0 + u_jet*u_jet );
      const double u_jet_perp = u_jet*Jet_Sine;
      const double u_jet_par  = u_jet*Jet_Cosine;

//    set fluid variable inside source
      PriReal[DENS        ] = (real)Jet_Density;
      PriReal[mom_axis    ] = (real)u_jet_par;
      PriReal[mom_x       ] = (real)u_jet_perp*cos_phi;
      PriReal[mom_y       ] = (real)u_jet_perp*sin_phi;
      PriReal[ENGY        ] = (real)Amb_Pressure;
      PriReal[JetFieldIdx ] = (real)Jet_Density;
      PriReal[ICMFieldIdx ] = (real)0.0;
      PriReal[LobeFieldIdx] = (real)0.0;
      PriReal[IntFieldIdx ] = (real)0.0;

      Hydro_Pri2Con( PriReal, BVal, false, PassiveNorm_NVar, PassiveNorm_VarIdx,
                     EoS_DensPres2Eint_CPUPtr, EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                     EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

   } // if ( Jet_Fire  &&  x <= 1.0 )

   else
   {

      for (int v=0; v<NVar_Flu; v++) {
         TFluVarIdx = TFluVarIdxList[v];

         if ( TFluVarIdx == mom_axis )
             BVal[TFluVarIdx] = MIN( Array3D[v][k][j_ref][i_ref], 0.0 );
         else
	     BVal[TFluVarIdx] = Array3D[v][k][j_ref][i_ref];
     }

   } // if ( Jet_Fire  &&  x <= 1.0 ) ... else ...

} // FUNCTION : JetBC



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_JetDeflect
// Description :  Add the problem-specific fields
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddField() for each of the problem-specific field:
//                   --> Field label sent to AddField() will be used as the output name of the field
//                   --> Field index returned by AddField() can be used to access the field data
//                3. Pre-declared field indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewField_JetDeflect()
{

  if ( JetFieldIdx == Idx_Undefined )
    JetFieldIdx = AddField( "JetField", FIXUP_FLUX_YES, FIXUP_REST_YES,
                            NORMALIZE_YES, INTERP_FRAC_NO );
  if ( ICMFieldIdx == Idx_Undefined )
    ICMFieldIdx = AddField( "ICMField", FIXUP_FLUX_YES, FIXUP_REST_YES,
                            NORMALIZE_YES, INTERP_FRAC_NO );
  if ( LobeFieldIdx == Idx_Undefined )
    LobeFieldIdx = AddField( "LobeField", FIXUP_FLUX_YES, FIXUP_REST_YES,
                             NORMALIZE_YES, INTERP_FRAC_NO );
  if ( IntFieldIdx == Idx_Undefined )
    IntFieldIdx = AddField( "IntField", FIXUP_FLUX_YES, FIXUP_REST_YES,
                            NORMALIZE_YES, INTERP_FRAC_NO );

} // FUNCTION : AddNewField_JetICMWall
#endif // #if ( MODEL == HYDRO )


//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_JetDeflect
// Description :  User-defined flag criteria
//
// Note        :  1. Invoked by Flag_Check() using the function pointer "Flag_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k     : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv        : Refinement level of the target patch
//                PID       : ID of the target patch
//                Threshold : User-provided threshold for the flag operation, which is loaded from the
//                            file "Input__Flag_User"
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_JetDeflect( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   const real (*Dens)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS];         // density
   const real (*JetD)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[JetFieldIdx];  // jet density
   const real (*ICMD)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ICMFieldIdx];  // icm density

   /*
   const double dh     = amr->dh[lv];
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
   const double dr[3]  = { Pos[0]-amr->BoxCenter[0], Pos[1]-amr->BoxCenter[1], Pos[2]-amr->BoxCenter[2] };
   const double Radius = sqrt( SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]) );
   */

   const real Xjet = JetD[k][j][i] / Dens[k][j][i];
   const real XICM = ICMD[k][j][i] / Dens[k][j][i];
      
   bool Flag = ( Xjet > Threshold[0] ) || ( XICM > Threshold[0] ) ;

   return Flag;

} // FUNCTION : Flag_JetDeflect


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_JetDeflect
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_JetDeflect()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Init_Field_User_Ptr      = AddNewField_JetDeflect;
   Flag_User_Ptr            = Flag_JetDeflect;
   Flag_Region_Ptr          = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   BC_User_Ptr              = JetBC;
   Flu_ResetByUser_Func_Ptr = NULL;
   Output_User_Ptr          = NULL;
   Aux_Record_User_Ptr      = NULL;
   End_User_Ptr             = NULL;
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_JetDeflect
