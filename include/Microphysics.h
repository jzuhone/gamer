#ifndef __MICROPHYSICS_H__
#define __MICROPHYSICS_H__

// include "Macro.h" and "Typedef.h" here since the header "GAMER.h" is NOT included in GPU solvers
#ifdef __CUDACC__
# include "Macro.h"
# include "Typedef.h"
#else
# include "GAMER.h"
#endif

#ifdef VISCOSITY

# define CONSTANT_VISCOSITY 1
# define SPITZER_VISCOSITY 2
# define VISCOSITY_KINETIC_COEFF 1
# define VISCOSITY_DYNAMIC_COEFF 2
# define ISOTROPIC_VISCOSITY 1
# define ANISOTROPIC_VISCOSITY 2

real FreqPrefactor;

#endif // #ifdef VISCOSITY

#endif // #ifndef __MICROPHYSICS_H__
