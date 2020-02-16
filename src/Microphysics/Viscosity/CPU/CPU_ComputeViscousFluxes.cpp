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
void Hydro_ComputeViscousFluxes( const real g_FC_Var [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                 const real Flux_1Face[NCOMP_TOTAL_PLUS_MAG],
                                 const int i_flux, const int j_flux, const int k_flux,
                                 const int i_fc, const int j_fc, const int k_fc, 
                                 const int d, const real dt, const real dh, const double Time )
{

    real vx, vy, vz;

    mu = Hydro_ComputeViscosity( fluid, Gamma_m1, MinPres );

#ifdef MHD

#else 

//  Isotropic viscosity

    switch ( d ) {
        case 0:
            break;
        case 1:
            break;
        case 2:
            break; 
    }

#endif // #ifdef MHD
   
#  ifdef __CUDACC__
    __syncthreads();
#  endif


} // FUNCTION : Hydro_ComputeViscousFluxes

#endif // #if ( ( MODEL == HYDRO ) && defined VISCOSITY )

#endif // #ifndef __CPU_VISCOUS_FLUXES__
