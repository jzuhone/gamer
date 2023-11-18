#ifndef __CUFLU_ADDCONDUCTIVEFLUX__
#define __CUFLU_ADDCONDUCTIVEFLUX__



#include "CUFLU.h"

#ifdef CONDUCTION



// external functions
#ifdef __CUDACC__

# include "CUFLU_ComputeConduction.cu"

#else // #ifdef __CUDACC__

void Hydro_ComputeConduction( real &cond_kappa, real &cond_chi, const MicroPhy_t *MicroPhy, 
                              const real Dens, const real Temp )

#endif // #ifdef __CUDACC__ ... else ...


// internal funcitons
GPU_DEVICE
static real MC_limiter( const real a, const real b );
GPU_DEVICE
static real minmod( const real a, const real b );




//-----------------------------------------------------------------------------------------
// Function    : Hydro_AddConductiveFlux_HalfStep
//
// Description : Compute the half-step conductive fluxes
//
// Note        : 1. Must enable CONDUCTION
//               2. Must enable MHD for anisotropic thermal conduction
//               3. Invoked by CPU/CUFLU_FluidSolver_MHM()
//
// Reference   : 
//
// Parameter   : g_Con_Var   : Array storing the input cell-centered temperature
//               g_Flux_Half : Array with hydrodynamic fluxes for adding the conductive fluxes
//               g_FC_B      : Array storing the input face-centered B field
//               g_CC_B      : Array storing the input cell-centered B field
//               dh          : Cell size
//               MicroPhy    : Microphysics object
//
// Return      : g_Flux_Half[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_AddConductiveFlux_HalfStep( const real Dens[], const real Temp[],
                                       real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                       const real g_FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                                       const real g_CC_B[][ CUBE(FLU_NXT) ],
                                       const real dh, const MicroPhy_t *MicroPhy )
{

   const int  didx_cvar[3] = { 1, FLU_NXT, SQR(FLU_NXT) };
   const int  flux_offset  = 1;  // skip the additional fluxes along the transverse directions for computing the CT electric field
   const real _dh          = (real)1.0 / dh;

   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

      int sizeB_i, sizeB_j, stride_fc_B;
      int size_i, size_j, size_k;
      int i_offset, j_offset, k_offset;

      switch ( d )
      {
         case 0 : size_i   = N_HF_FLUX-1;              size_j   = N_HF_FLUX-2*flux_offset;  size_k      = N_HF_FLUX-2*flux_offset;
                  i_offset = 0;                        j_offset = flux_offset;              k_offset    = flux_offset;
                  sizeB_i  = FLU_NXT_P1;               sizeB_j  = FLU_NXT;                  stride_fc_B = 1;
                  break;

         case 1 : size_i   = N_HF_FLUX-2*flux_offset;  size_j   = N_HF_FLUX-1;              size_k      = N_HF_FLUX-2*flux_offset;
                  i_offset = flux_offset;              j_offset = 0;                        k_offset    = flux_offset;
                  sizeB_i  = FLU_NXT;                  sizeB_j  = FLU_NXT_P1;               stride_fc_B = FLU_NXT;
                  break;

         case 2 : size_i   = N_HF_FLUX-2*flux_offset;  size_j   = N_HF_FLUX-2*flux_offset;  size_k      = N_HF_FLUX-1;
                  i_offset = flux_offset;              j_offset = flux_offset;              k_offset    = 0;
                  sizeB_i  = FLU_NXT;                  sizeB_j  = FLU_NXT;                  stride_fc_B = SQR(FLU_NXT);
                  break;
      } // switch ( d )

      const int size_ij = size_i*size_j;

      CGPU_LOOP( idx, size_i*size_j*size_k )
      {
//       flux index
         const int i_flux   = idx % size_i           + i_offset;
         const int j_flux   = idx % size_ij / size_i + j_offset;
         const int k_flux   = idx / size_ij          + k_offset;
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, N_HF_FLUX, N_HF_FLUX );

//       conserved variable and cell-centered magnetic field index
         const int i_cvar   = i_flux;
         const int j_cvar   = j_flux;
         const int k_cvar   = k_flux;
         const int idx_cvar = IDX321( i_cvar, j_cvar, k_cvar, FLU_NXT, FLU_NXT );

//       face-centered magnetic field index
         const int idx_fc_B = IDX321( i_cvar, j_cvar, k_cvar, sizeB_i, sizeB_j ) + stride_fc_B;


//       1. get the conductivity
//           --> for non-constant conduction coefficients, take the spatial average along the normal direction
//               to get the face-centered coefficients
         real kappa_l, kappa_r, kappa, chi;
         Hydro_ComputeConduction( kappa_l, cond_chi, MicroPhy, Dens[ idx_cvar ], 
                                  Temp[ idx_cvar ] );
         Hydro_ComputeConduction( kappa_r, cond_chi, MicroPhy, Dens[ idx_cvar + didx_cvar[d] ], 
                                  Temp[ idx_cvar + didx_cvar[d] ] );
         kappa = 0.5*(kappa_l + kappa_r);

//       2. compute the mean magnetic field
//       ---------------------
//       |         |         |
//       |    ^    |    ^    |
//       -----1---------2-----
//       |    |    |    |    |
//       |         |         |
//       |  i j k -->        |
//       |         |         |
//       |    ^    |    ^    |
//       -----3---------4-----
//       |    |    |    |    |
//       |         |         |
//       ---------------------
         real B_N_mean, B_T1_mean, B_T2_mean, B_amp;
         B_N_mean  =             g_FC_B[    d][ idx_fc_B                ];
         B_T1_mean = (real)0.5*( g_CC_B[TDir1][ idx_cvar                ] +
                                 g_CC_B[TDir1][ idx_cvar + didx_cvar[d] ]   );
         B_T2_mean = (real)0.5*( g_CC_B[TDir2][ idx_cvar                ] +
                                 g_CC_B[TDir2][ idx_cvar + didx_cvar[d] ]   );
         B_amp     = SQRT( SQR(B_N_mean) + SQR(B_T1_mean) + SQR(B_T2_mean) );

//       normalize magnetic field
         B_N_mean  /= B_amp;
         B_T1_mean /= B_amp;
         B_T2_mean /= B_amp;


//       3. compute temperature slope
//       ---------------------
//       |         |         |
//       ----bl--------br-----
//       |         |         |
//       |      N_slope      |
//       |         |         |
//       ----al--------ar-----
//       |         |         |
//       ---------------------
         real N_slope, T1_slope, T2_slope;
         real al, bl, ar, br;

//       normal direction
         N_slope = ( Temp[ idx_cvar + didx_cvar[d] ] - Temp[ idx_cvar ] ) * _dh;

//       transverse direction 1
         al = Temp[ idx_cvar                                   ] -
              Temp[ idx_cvar                - didx_cvar[TDir1] ];
         bl = Temp[ idx_cvar                + didx_cvar[TDir1] ] -
              Temp[ idx_cvar                                   ];
         ar = Temp[ idx_cvar + didx_cvar[d]                    ] -
              Temp[ idx_cvar + didx_cvar[d] - didx_cvar[TDir1] ];
         br = Temp[ idx_cvar + didx_cvar[d] + didx_cvar[TDir1] ] -
              Temp[ idx_cvar + didx_cvar[d]                    ];
         T1_slope = (  MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) )  ) * _dh;

//       transverse direction 2
         al = Temp[ idx_cvar                                   ] -
              Temp[ idx_cvar                - didx_cvar[TDir2] ];
         bl = Temp[ idx_cvar                + didx_cvar[TDir2] ] -
              Temp[ idx_cvar                                   ];
         ar = Temp[ idx_cvar + didx_cvar[d]                    ] -
              Temp[ idx_cvar + didx_cvar[d] - didx_cvar[TDir2] ];
         br = Temp[ idx_cvar + didx_cvar[d] + didx_cvar[TDir2] ] -
              Temp[ idx_cvar + didx_cvar[d]                    ];
         T2_slope = (  MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) )  ) * _dh;


//       4. compute conductive flux
         real Flux_Total, gradient;

         gradient   = -B_N_mean*( B_N_mean*N_slope + B_T1_mean*T1_slope + B_T2_mean*T2_slope );
         Flux_Total = kappa*gradient;

//       6. flux add-up
         g_Flux_Half[d][ENGY][idx_flux] += Flux_Total;

      } // CGPU_LOOP( idx, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_AddConductiveFlux_HalfStep



//-----------------------------------------------------------------------------------------
// Function    : Hydro_AddConductiveFlux_FullStep
//
// Description : Compute the full-step conductive fluxes
//
// Note        : 1. Must enable CONDUCTION
//               2. Must enable MHD for anisotropic thermal conduction
//               3. Invoked by CPU/CUFLU_FluidSolver_MHM()
//
// Reference   : Yang et al., ApJ 761, 185 (2012); doi:10.1088/0004-637X/761/2/185
//
// Parameter   : g_PriVar_Half : Array storing the input cell-centered, half-step primitive fluid variables
//               g_FC_Flux     : Array with hydrodynamic fluxes for adding the conductive fluxes
//               g_FC_B_Half   : Array storing the input face-centered, half-step magnetic field
//               NFlux         : Stride for accessing g_FC_Flux[]
//               dh            : Cell size
//               MicroPhy      : Microphysics object
//
// Return      : g_FC_Flux[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_AddConductiveFlux_FullStep( const real Dens[], const real Temp[],
                                       real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                       const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                       const int NFlux, const real dh, const MicroPhy_t *MicroPhy )
{

   const int  didx_half[3] = { 1, N_HF_VAR, SQR(N_HF_VAR) };
   const int  mag_offset   = ( N_HF_VAR - PS2 ) / 2;
   const int  half_offset  = ( N_HF_VAR - NFlux ) / 2;
   const int  flux_offset  = 1;  // skip the additional fluxes along the transverse directions for computing the CT electric field
   const real _dh          = (real)1.0 / dh;

   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

      int sizeB_i, sizeB_j, sizeB_k;
      int mag_offset_i, mag_offset_j, mag_offset_k;
      int size_i, size_j, size_k;
      int i_offset, j_offset, k_offset;

      switch ( d )
      {
         case 0 : size_i       = NFlux-1;               size_j        = NFlux-2*flux_offset;  size_k        = NFlux-2*flux_offset;
                  sizeB_i      = N_HF_VAR+1;            sizeB_j       = N_HF_VAR;             sizeB_k       = N_HF_VAR;
                  mag_offset_i = mag_offset;            mag_offset_j  = mag_offset-1;         mag_offset_k  = mag_offset-1;
                  i_offset     = 0;                     j_offset      = flux_offset;          k_offset      = flux_offset;
                  break;

         case 1 : size_i        = NFlux-2*flux_offset;  size_j        = NFlux-1;              size_k        = NFlux-2*flux_offset;
                  sizeB_i       = N_HF_VAR;             sizeB_j       = N_HF_VAR+1;           sizeB_k       = N_HF_VAR;
                  mag_offset_i  = mag_offset-1;         mag_offset_j  = mag_offset;           mag_offset_k  = mag_offset-1;
                  i_offset      = flux_offset;          j_offset      = 0;                    k_offset      = flux_offset;
                  break;

         case 2 : size_i        = NFlux-2*flux_offset;  size_j        = NFlux-2*flux_offset;  size_k        = NFlux-1;
                  sizeB_i       = N_HF_VAR;             sizeB_j       = N_HF_VAR;             sizeB_k       = N_HF_VAR+1;
                  mag_offset_i  = mag_offset-1;         mag_offset_j  = mag_offset-1;         mag_offset_k  = mag_offset;
                  i_offset      = flux_offset;          j_offset      = flux_offset;          k_offset      = 0;
                  break;
      } // switch ( d )

      const int stride_fc_BT1[3] = { 1, sizeB_k, sizeB_k*sizeB_i };
      const int stride_fc_BT2[3] = { 1, sizeB_j, sizeB_j*sizeB_k };
      const int size_ij          = size_i*size_j;

      CGPU_LOOP( idx, size_i*size_j*size_k )
      {
//       flux index
         const int i_flux     = idx % size_i           + i_offset;
         const int j_flux     = idx % size_ij / size_i + j_offset;
         const int k_flux     = idx / size_ij          + k_offset;
         const int idx_flux   = IDX321( i_flux, j_flux, k_flux, NFlux, NFlux );

//       half-step primitive variable index
         const int i_half     = i_flux + half_offset;
         const int j_half     = j_flux + half_offset;
         const int k_half     = k_flux + half_offset;
         const int idx_half   = IDX321( i_half, j_half, k_half, N_HF_VAR, N_HF_VAR );

//       magnetic field indices
         const int i_fc       = i_flux + mag_offset_i;
         const int j_fc       = j_flux + mag_offset_j;
         const int k_fc       = k_flux + mag_offset_k;
         const int idx_fc_BN  = IDX321( i_fc, j_fc, k_fc, sizeB_i, sizeB_j );
         const int idx_fc_BT1 = IDX321( i_fc, j_fc, k_fc, sizeB_k, sizeB_i );
         const int idx_fc_BT2 = IDX321( i_fc, j_fc, k_fc, sizeB_j, sizeB_k );


//       1. get the conductivity
//           --> for non-constant conduction coefficients, take the spatial average along the normal direction
//               to get the face-centered coefficients
         real kappa_l, kappa_r, kappa, chi;
         Hydro_ComputeConduction( kappa_l, chi, MicroPhy, Dens[ idx_half ], 
                                  Temp[ idx_half ] );
         Hydro_ComputeConduction( kappa_r, chi, MicroPhy, Dens[ idx_half + didx_half[d] ], 
                                  Temp[ idx_half + didx_half[d] ] );
         kappa = 0.5*(kappa_l + kappa_r);

//       2. compute the mean magnetic field
//       ---------------------
//       |         |         |
//       |    ^    |    ^    |
//       -----1---------2-----
//       |    |    |    |    |
//       |         |         |
//       |  i j k -->        |
//       |         |         |
//       |    ^    |    ^    |
//       -----3---------4-----
//       |    |    |    |    |
//       |         |         |
//       ---------------------
         real B_N_mean, B_T1_mean, B_T2_mean, B_amp;
         B_N_mean  =              g_FC_B_Half[    d][ idx_fc_BN                                            ];
         B_T1_mean = (real)0.25*( g_FC_B_Half[TDir1][ idx_fc_BT1                                           ] +
                                  g_FC_B_Half[TDir1][ idx_fc_BT1                    + stride_fc_BT1[TDir1] ] +
                                  g_FC_B_Half[TDir1][ idx_fc_BT1 - stride_fc_BT1[d]                        ] +
                                  g_FC_B_Half[TDir1][ idx_fc_BT1 - stride_fc_BT1[d] + stride_fc_BT1[TDir1] ]   );
         B_T2_mean = (real)0.25*( g_FC_B_Half[TDir2][ idx_fc_BT2                                           ] +
                                  g_FC_B_Half[TDir2][ idx_fc_BT2                    + stride_fc_BT2[TDir2] ] +
                                  g_FC_B_Half[TDir2][ idx_fc_BT2 - stride_fc_BT2[d]                        ] +
                                  g_FC_B_Half[TDir2][ idx_fc_BT2 - stride_fc_BT2[d] + stride_fc_BT2[TDir2] ]   );
         B_amp     = SQRT( SQR(B_N_mean) + SQR(B_T1_mean) + SQR(B_T2_mean) );

//       normalize magnetic field
         B_N_mean  /= B_amp;
         B_T1_mean /= B_amp;
         B_T2_mean /= B_amp;


//       3. compute temperature slope
//       ---------------------
//       |         |         |
//       ----bl--------br-----
//       |         |         |
//       |      N_slope      |
//       |         |         |
//       ----al--------ar-----
//       |         |         |
//       ---------------------
         real N_slope, T1_slope, T2_slope;
         real al, bl, ar, br;

//       normal direction
         N_slope = ( Temp[ idx_half + didx_half[d] ] - Temp[ idx_half ] ) * _dh;

//       transverse direction 1
         al = Temp[ idx_half                                   ] -
              Temp[ idx_half                - didx_half[TDir1] ];
         bl = Temp[ idx_half                + didx_half[TDir1] ] -
              Temp[ idx_half                                   ];
         ar = Temp[ idx_half + didx_half[d]                    ] -
              Temp[ idx_half + didx_half[d] - didx_half[TDir1] ];
         br = Temp[ idx_half + didx_half[d] + didx_half[TDir1] ] -
              Temp[ idx_half + didx_half[d]                    ];
         T1_slope = (  MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) )  ) * _dh;

//       transverse direction 2
         al = Temp[ idx_half                                   ] -
              Temp[ idx_half                - didx_half[TDir2] ];
         bl = Temp[ idx_half                + didx_half[TDir2] ] -
              Temp[ idx_half                                   ];
         ar = Temp[ idx_half + didx_half[d]                    ] -
              Temp[ idx_half + didx_half[d] - didx_half[TDir2] ];
         br = Temp[ idx_half + didx_half[d] + didx_half[TDir2] ] -
              Temp[ idx_half + didx_half[d]                    ];
         T2_slope = (  MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) )  ) * _dh;


//       4. compute CR diffusive flux
         real Flux_Total, gradient;

         gradient   = -B_N_mean*( B_N_mean*N_slope + B_T1_mean*T1_slope + B_T2_mean*T2_slope );
         Flux_Total = kappa*gradient;


//       5. flux add-up
         g_FC_Flux[d][ENGY][idx_flux] += Flux_Total;

      } // CGPU_LOOP( idx, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_AddConductiveFlux_FullStep



//-----------------------------------------------------------------------------------------
// Function    : MC_limiter
// Description : Monotonized central (MC) slope limiter
//
// Parameter   : a, b : Input slopes
//
// Return      : Limited slope
//-----------------------------------------------------------------------------------------
GPU_DEVICE
static real MC_limiter( const real a, const real b )
{

   return minmod( (real)2.0*minmod(a,b), (real)0.5*(a+b) );

} // FUNCTION : MC_limiter



//-----------------------------------------------------------------------------------------
// Function    : minmod
// Description : Minmod slope limiter
//
// Parameter   : a, b : Input slopes
//
// Return      : Limited slope
//-----------------------------------------------------------------------------------------
GPU_DEVICE
static real minmod( const real a, const real b )
{

   if      ( a > (real)0.0  &&  b > (real)0.0 )    return FMIN(a, b);
   else if ( a < (real)0.0  &&  b < (real)0.0 )    return FMAX(a, b);
   else                                            return (real)0.0;

} // FUNCTION : minmod



#endif // #ifdef CONDUCTION



#endif // #ifndef __CUFLU_ADDCONDUCTIVEFLUX__