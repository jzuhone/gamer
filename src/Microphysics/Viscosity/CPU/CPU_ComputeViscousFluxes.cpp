#ifndef __CPU_VISCOUS_FLUXES__
#define __CPU_VISCOUS_FLUXES__

#include "Microphysics.h"

#if ( ( MODEL == HYDRO ) && defined VISCOSITY )

// external functions
#ifdef __CUDACC__
# include "CUFLU_Viscosity.cu"
#else // #ifdef __CUDACC__

#endif // #ifdef __CUDACC__ ... else ...

GPU_DEVICE
void Hydro_ComputeIsoViscousFluxes( const real g_FC_Var [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                                    real g_FC_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                                    const real dt, const real dh, const double Time )
{


   
#  ifdef __CUDACC__
    __syncthreads();
#  endif

} // FUNCTION : Hydro_ComputeIsoViscousFluxes

#ifdef MHD

GPU_DEVICE
void Hydro_ComputeAnisoViscousFluxes( const real g_FC_Var [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                                      real g_FC_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                                      const real dt, const real dh, const double Time )
{


   
#  ifdef __CUDACC__
    __syncthreads();
#  endif

} // FUNCTION : Hydro_ComputeAnisoViscousFluxes

#endif // #ifdef MHD

#endif // #if ( ( MODEL == HYDRO ) && defined VISCOSITY )

#endif // #ifndef __CPU_VISCOUS_FLUXES__
