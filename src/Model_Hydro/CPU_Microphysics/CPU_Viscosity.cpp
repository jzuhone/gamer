#ifndef __CUFLU_VISCOSITY__
#define __CUFLU_VISCOSITY__

#include "CUFLU.h"

#if ( ( MODEL == HYDRO  || MODEL == MHD ) &&  defined VISCOSITY )

// external functions
#ifdef __CUDACC__
# include "CUFLU_Viscosity.cu"
#else // #ifdef __CUDACC__

#endif // #ifdef __CUDACC__ ... else ...

GPU_DEVICE
void Hydro_ViscosityInit ( )
{
    if ( VISCOSITY_TYPE == SPITZER_VISCOSITY ) {
        real lnLambda = 30; // Coulomb logarithm for Spitzer viscosity
        // This prefactor is in CGS and we must convert it to code units
        // To avoid precision errors, we do this one step at a time
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

} // FUNCTION : Hydro_ViscosityInit

GPU_DEVICE
void Hydro_ComputeViscosity( real nu[ CUBE(N_FC_VAR) ],
                             const real Flu_Array[NCOMP_FLUID][ CUBE(PS1) ],
                             const real Gamma_m1, const real MinPres )
{
    const bool CheckMinPres_Yes = true;

    CGPU_LOOP( t, CUBE(PS1) )
    {
        real fluid[NCOMP_FLUID], _Rho;

        for (int v=0; v<NCOMP_FLUID; v++)   fluid[v] = g_Flu_Array[v][t];

        _Rho  = (real)1.0 / fluid[DENS];

        if ( VISCOSITY_TYPE == CONSTANT_VISCOSITY ) {

            // Constant viscosity
            if ( VISCOSITY_COEFF_TYPE == VISCOSITY_KINETIC_COEFF ) {

                nu[t] = (real)VISCOSITY_COEFF;

            } else if ( VISCOSITY_COEFF_TYPE == VISCOSITY_DYNAMIC_COEFF ) {

                nu[t] = (real)VISCOSITY_COEFF*_Rho;

            }

        } else if ( VISCOSITY_TYPE == SPITZER_VISCOSITY ) { 

            // Spitzer viscosity
            real Pres, Temp, Freq_ii;

            Pres = Hydro_GetPressure( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], 
                                      fluid[ENGY], Gamma_m1, CheckMinPres_Yes, MinPres );
            
            Temp = Pres*_Rho;

            Freq_ii = FreqPrefactor*fluid[DENS]*POW( Temp, (real)-1.5 );

            nu[t] = VISCOSITY_SPITZER_FRACTION*0.96*Pres/Freq_ii;

        }

        nu[t] = FMIN( FMAX( nu[t], VISCOSITY_COEFF_MIN ), VISCOSITY_COEFF_MAX );

    } // CGPU_LOOP( t, CUBE(PS1) )

} // FUNCTION : Hydro_ComputeViscosity

GPU_DEVICE
void Hydro_ComputeIsoViscousFluxes( const real g_FC_Var [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                                    real g_FC_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                                    const real dt, const real dh, const double Time )
{


   
#  ifdef __CUDACC__
    __syncthreads();
#  endif

} // FUNCTION : Hydro_ComputeIsoViscousFluxes

#if ( MODEL == MHD )

GPU_DEVICE
void Hydro_ComputeAnisoViscousFluxes( const real g_FC_Var [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                                      real g_FC_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                                      const real dt, const real dh, const double Time )
{


   
#  ifdef __CUDACC__
    __syncthreads();
#  endif

} // FUNCTION : Hydro_ComputeAnisoViscousFluxes

#endif // #if ( MODEL == MHD )

//-----------------------------------------------------------------------------------------
// Function    :  CPU/CUFLU_dtSolver_HydroCFL
// Description :  Estimate the evolution time-step (dt) from the CFL condition of viscosity
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. time-step is estimated by the stability criterion from the von Neumann stability analysis
//                   --> CFL condition
//                3. Arrays with a prefix "g_" are stored in the global memory of GPU
//
// Parameter   :  g_dt_Array  : Array to store the minimum dt in each target patch
//                g_Flu_Array : Array storing the prepared fluid data of each target patch
//                NPG         : Number of target patch groups (for CPU only)
//                dh          : Cell size
//                Safety      : dt safety factor
//                Gamma       : Ratio of specific heats
//                MinPres     : Minimum allowed pressure
//
// Return      :  g_dt_Array
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUFLU_dtSolver_ViscCFL( real g_dt_Array[], const real g_Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ],
                              const real dh, const real Safety, const real Gamma, const real MinPres )
#else
void CPU_dtSolver_ViscCFL  ( real g_dt_Array[], const real g_Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ], const int NPG,
                              const real dh, const real Safety, const real Gamma, const real MinPres )
#endif
{

   const real Gamma_m1         = Gamma - (real)1.0;
   const real dh2Safety         = Safety*0.5*dh*dh;

// loop over all patches
// --> CPU/GPU solver: use different (OpenMP threads) / (CUDA thread blocks)
//                     to work on different patches
#  ifdef __CUDACC__
   const int p = blockIdx.x;
#  else
#  pragma omp parallel for schedule( runtime )
   for (int p=0; p<8*NPG; p++)
#  endif
   {
      real MaxCFL=(real)0.0;

      Hydro_ComputeViscosity( nu, g_Flu_Array[p], Gamma_m1, MinPres ); 

      CGPU_LOOP( t, CUBE(PS1) )
      {

         MaxCFL = FMAX( nu[t], MaxCFL );

      } // CGPU_LOOP( t, CUBE(PS1) )

//    perform parallel reduction to get the maximum CFL speed in each thread block
//    --> store in the thread 0
#     ifdef __CUDACC__
#     ifdef DT_FLU_USE_SHUFFLE
      MaxCFL = BlockReduction_Shuffle ( MaxCFL );
#     else
      MaxCFL = BlockReduction_WarpSync( MaxCFL );
#     endif
      if ( threadIdx.x == 0 )
#     endif // #ifdef __CUDACC__
      g_dt_Array[p] = dh2Safety/MaxCFL;

   } // for (int p=0; p<8*NPG; p++)

} // FUNCTION : CPU/CUFLU_dtSolver_ViscCFL


#endif // #if ( ( MODEL == HYDRO || MODEL == MHD )  &&  defined VISCOSITY )

#endif // #ifndef __CUFLU_VISCOSITY__

