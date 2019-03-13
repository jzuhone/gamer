#ifndef __CUFLU_VISCOSITY__
#define __CUFLU_VISCOSITY__

#include "Microphysics.h"

#if ( ( MODEL == HYDRO  || MODEL == MHD ) &&  defined VISCOSITY )

// external functions
#ifdef __CUDACC__
# include "CUFLU_Viscosity.cu"
#else // #ifdef __CUDACC__

#endif // #ifdef __CUDACC__ ... else ...

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

#endif // #if ( ( MODEL == HYDRO || MODEL == MHD )  &&  defined VISCOSITY )

#endif // #ifndef __CUFLU_VISCOSITY__

