#ifndef __COMPUTE_VISCOSITY__
#define __COMPUTE_VISCOSITY__

#include "Microphysics.h"

#if ( ( MODEL == HYDRO ) && defined VISCOSITY )

GPU_DEVICE
real Hydro_ComputeViscosity( const real fluid[NCOMP_FLUID],
                             const real Gamma_m1, const real MinPres )
{
    const bool CheckMinPres_Yes = true;

    real nu, _Rho;

    _Rho  = (real)1.0 / fluid[DENS];

    if ( VISCOSITY_TYPE == CONSTANT_VISCOSITY ) {

        // Constant viscosity
        if ( VISCOSITY_COEFF_TYPE == VISCOSITY_KINETIC_COEFF ) {

            nu = (real)VISCOSITY_COEFF;

        } else if ( VISCOSITY_COEFF_TYPE == VISCOSITY_DYNAMIC_COEFF ) {

            nu = (real)VISCOSITY_COEFF*_Rho;

        }

    } else if ( VISCOSITY_TYPE == SPITZER_VISCOSITY ) { 

        // Spitzer viscosity
        real Pres, Temp, Freq_ii;

        Pres = Hydro_GetPressure( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], 
                                  fluid[ENGY], Gamma_m1, CheckMinPres_Yes, MinPres );
        
        Temp = Pres*_Rho;

        Freq_ii = FreqPrefactor*fluid[DENS]*POW( Temp, (real)-1.5 );

        nu = VISCOSITY_SPITZER_FRACTION*0.96*Pres/Freq_ii;

    }

    nu = FMIN( FMAX( nu, VISCOSITY_COEFF_MIN ), VISCOSITY_COEFF_MAX );

    return nu;

} // FUNCTION : Hydro_ComputeViscosity

#endif // #if ( ( MODEL == HYDRO ) && defined VISCOSITY )

#endif // #ifndef __COMPUTE_VISCOSITY__

