#include "GAMER.h"

#ifdef SUPPORT_HDF5
#include "hdf5.h"
#endif

static int nAx, nAy, nAz;
static double Axmin, Aymin, Azmin;
static double Adx, Ady, Adz;
static double *Axcoord, *Aycoord, *Azcoord;

// Grid metadata (coordinate arrays, spacing, lower edge) needed by VecPot_Interp.
// A GridMeta can describe either the root-level coarse grid (Axcoord/Adx/Axmin, etc.)
// or a single refinement patch (BFieldPatch, below) -- VecPot_Interp itself doesn't
// care which.
struct GridMeta
{
   const double *xcoord, *ycoord, *zcoord;
   double        dx, dy, dz;
   double        xmin, ymin, zmin;
};

// A single locally-refined B-field vector-potential patch, optionally present in the
// input file as flat HDF5 groups "patch_0000", "patch_0001", ... (see NumBPatches
// below). Optional: files with no such groups (a plain uniform-mesh vector-potential
// file, with no refinement patches at all) behave exactly as before -- see NumBPatches
// below.
//
// Every patch -- however deeply nested -- is stored flat in BPatches, one entry per
// "patch_NNNN" group, carrying no parent/child reference of any kind: just its own box
// (coordinates and field data -- already the final, fully-blended vector potential,
// not a correction needing further combination with anything coarser) plus "level",
// its own nesting depth (1 for a patch nested directly inside the coarse grid, 2 for
// one nested inside a level-1 patch, etc.). Since every patch's data already stands on
// its own, evaluating a point just means picking whichever patch *at the finest level*
// contains it and reading that patch's own value directly -- see VecPot_EvalComponent
// -- with no blending step and no need to know which coarser patch a given one sits
// inside. This keeps the on-disk format itself down to "a flat list of boxes, each
// tagged with a level and carrying its own self-contained field data", producible by
// any tool that can write such boxes out, not tied to any particular IC generator.
struct BFieldPatch
{
   int         nx, ny, nz;
   double      dx, dy, dz;
   double      xmin, ymin, zmin;
   double     *xcoord, *ycoord, *zcoord;
   double     *Ax, *Ay, *Az;
   int         level;
};

static int         NumBPatches     = 0;
static int         MaxBPatchLevel  = 0;
static BFieldPatch *BPatches       = NULL;

double TSC_Weight( const double x );
double VecPot_Interp( const double field[], const double xx, const double yy,
                      const double zz, const int fdims[], const int fbegin[],
                      const GridMeta &grid );
bool   VecPot_PointInPatch( const BFieldPatch &patch, const double xx, const double yy, const double zz );
double VecPot_EvalComponent( const double coarse_field[], const int fdims[], const int fbegin[],
                             const double xx, const double yy, const double zz, const int comp );
void   VecPot_FreePatch( BFieldPatch &patch );
#ifdef SUPPORT_HDF5
void VecPot_ReadField( hid_t mag_file_id, const int ibegin, const int jbegin,
                       const int kbegin, const int iend, const int jend,
                       const int kend, double Ax[], double Ay[], double Az[] );
void VecPot_ReadPatches( hid_t mag_file_id );
void VecPot_ReadPatchGroup( hid_t group_id, BFieldPatch &patch );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_Init_BField_ByVecPot_File
// Description :  Use the input uniform-mesh array (plus potential finer meshes) stored in
//                the file "B_Filename" to assign a magnetic vector potential to all real
//                patches on level "B_lv" and take the curl of this potential to compute
//                the magnetic field
//
// Note        :
//
// Parameter   :  B_lv : Target AMR level
//
// Return      :  amr->patch->magnetic
//-------------------------------------------------------------------------------------------------------
void MHD_Init_BField_ByVecPot_File( const int B_lv )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading the magnetic field from the input file ...\n" );


#  ifndef MHD
   Aux_Error( ERROR_INFO, "MHD must be enabled !!\n" );
#  endif

#  ifndef SUPPORT_HDF5
   Aux_Error( ERROR_INFO, "SUPPORT_HDF5 must be set to load a vector potential from a file !!\n" );
#  endif


   const char B_Filename[] = "B_IC";

   const double dh         = amr->dh[B_lv];

   double *Axf, *Ayf, *Azf;

   if ( !Aux_CheckFileExist(B_Filename) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", B_Filename );

// Open the magnetic field file and determine the dimensionality of the vector
// potential grid

#  ifdef SUPPORT_HDF5

   herr_t status;
   hid_t dataset, dataspace;
   hsize_t dims[3], maxdims[3];
   int ndim;

   hid_t mag_file_id = H5Fopen(B_Filename, H5F_ACC_RDONLY, H5P_DEFAULT);

   if ( B_lv == 0 ) {

     dataset = H5Dopen(mag_file_id, "magnetic_vector_potential_z", H5P_DEFAULT);

     dataspace = H5Dget_space(dataset);

     ndim = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

     if ( ndim != 3 ) Aux_Error( ERROR_INFO, "Incorrect dimensionality of vector potential ndim=%d !!\n", ndim );

     H5Sclose(dataspace);
     H5Dclose(dataset);

//   NOTE: Magnetic vector potential arrays are stored in column-major order,
//   i.e., Ax[nAx][nAy][nAz]

     nAx = dims[0];
     nAy = dims[1];
     nAz = dims[2];

//   Optionally read any locally-refined vector-potential patches (see BFieldPatch).
//   A plain uniform-mesh file simply has no "num_patches" attribute -- NumBPatches
//   stays 0 and everything below behaves exactly as it did before this feature
//   existed.
     VecPot_ReadPatches( mag_file_id );

     if ( MPI_Rank == 0 && NumBPatches > 0 )
        Aux_Message( stdout, "   Loading %d B-field refinement patch(es) from \"%s\" ...\n", NumBPatches, B_Filename );

   }

#  endif

// Read the coordinate information from the vector potential grid
   Axcoord = new double [ nAx ];
   Aycoord = new double [ nAy ];
   Azcoord = new double [ nAz ];

#  ifdef SUPPORT_HDF5

   dataset = H5Dopen(mag_file_id, "x", H5P_DEFAULT);
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, Axcoord);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load x-coordinate !!\n" );

   H5Dclose(dataset);

   dataset = H5Dopen(mag_file_id, "y", H5P_DEFAULT);
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, Aycoord);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load y-coordinate !!\n" );
   H5Dclose(dataset);

   dataset = H5Dopen(mag_file_id, "z", H5P_DEFAULT);
   status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, Azcoord);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load z-coordinate !!\n" );
   H5Dclose(dataset);

#  endif

// Cell spacing and left edge of vector potential grid

   Adx = Axcoord[1]-Axcoord[0];
   Ady = Aycoord[1]-Aycoord[0];
   Adz = Azcoord[1]-Azcoord[0];

   Axmin = Axcoord[0]-0.5*Adx;
   Aymin = Aycoord[0]-0.5*Ady;
   Azmin = Azcoord[0]-0.5*Adz;

   double Axmax = Axcoord[nAx-1]+0.5*Adx;
   double Aymax = Aycoord[nAy-1]+0.5*Ady;
   double Azmax = Azcoord[nAz-1]+0.5*Adz;

   if ( amr->BoxEdgeL[0] < Axmin+2*Adx || amr->BoxEdgeR[0] >= Axmax-2*Adx ||
        amr->BoxEdgeL[1] < Aymin+2*Ady || amr->BoxEdgeR[1] >= Aymax-2*Ady ||
        amr->BoxEdgeL[2] < Azmin+2*Adz || amr->BoxEdgeR[2] >= Azmax-2*Adz )
      Aux_Error( ERROR_INFO, "Input grid is smaller than the simulation domain !!" );

   double *Ax = new double [ CUBE(PS1+1) ];
   double *Ay = new double [ CUBE(PS1+1) ];
   double *Az = new double [ CUBE(PS1+1) ];

   double sample_res = pow(2, MAX_LEVEL-B_lv);
   double sample_fact = 1.0/((double)sample_res);

   for (int PID=0; PID<amr->NPatchComma[B_lv][1]; PID++) {

      double EdgeL[3];
      double EdgeR[3];

      for (int i=0; i<3; i++) {
         EdgeL[i] = amr->patch[0][B_lv][PID]->EdgeL[i];
         EdgeR[i] = amr->patch[0][B_lv][PID]->EdgeR[i];
      }

//    Compute the beginning and ending indices on the vector potential grid
//    +/- 1 are necessary because we will be computing derivatives

      int ibegin = (int)((EdgeL[0]-Axmin)/Adx)-2;
      int jbegin = (int)((EdgeL[1]-Aymin)/Ady)-2;
      int kbegin = (int)((EdgeL[2]-Azmin)/Adz)-2;

      int iend   = (int)((EdgeR[0]-Axmin)/Adx)+2;
      int jend   = (int)((EdgeR[1]-Aymin)/Ady)+2;
      int kend   = (int)((EdgeR[2]-Azmin)/Adz)+2;

      int nlocx  = iend-ibegin+1;
      int nlocy  = jend-jbegin+1;
      int nlocz  = kend-kbegin+1;

      int fdims[3] = { nlocx, nlocy, nlocz };
      int fbegin[3] = { ibegin, jbegin, kbegin };
      int nloc = nlocx*nlocy*nlocz;

//    Allocate for the data on the vector potential grid local to this patch and
//    read it from the file

      Axf = new double [nloc];
      Ayf = new double [nloc];
      Azf = new double [nloc];

#     ifdef SUPPORT_HDF5
      VecPot_ReadField( mag_file_id, ibegin, jbegin, kbegin,
                        iend, jend, kend, Axf, Ayf, Azf );
#     endif

//    Loop over the indices in this patch and interpolate the vector potential
//    to the current refinement level's resolution

      for (int k=0; k<PS1+1; k++) {  const double z0 = EdgeL[2] + k*dh;
      for (int j=0; j<PS1+1; j++) {  const double y0 = EdgeL[1] + j*dh;
      for (int i=0; i<PS1+1; i++) {  const double x0 = EdgeL[0] + i*dh;

         int idx = IDX321( i, j, k, PS1+1, PS1+1 );

         Ax[idx] = 0.0;
         Ay[idx] = 0.0;
         Az[idx] = 0.0;

         if ( i != PS1 ) {

            for ( int ii=0; ii<sample_res; ii++ ) {
               const double x = x0 + (ii+0.5)*dh*sample_fact;
               Ax[idx] += VecPot_EvalComponent( Axf, fdims, fbegin, x, y0, z0, 0 );
            }

         }

         if ( j != PS1 ) {

            for ( int jj=0; jj<sample_res; jj++ ) {
               const double y = y0 + (jj+0.5)*dh*sample_fact;
               Ay[idx] += VecPot_EvalComponent( Ayf, fdims, fbegin, x0, y, z0, 1 );
            }

         }

         if ( k != PS1 ) {

            for ( int kk=0; kk<sample_res; kk++ ) {
               const double z = z0 + (kk+0.5)*dh*sample_fact;
               Az[idx] += VecPot_EvalComponent( Azf, fdims, fbegin, x0, y0, z, 2 );
            }

         }

         Ax[idx] *= sample_fact;
         Ay[idx] *= sample_fact;
         Az[idx] *= sample_fact;

      }}}

//    Calculate Bx from vector potential
      for (int k=0; k<PS1;   k++) {
      for (int j=0; j<PS1;   j++) {
      for (int i=0; i<PS1+1; i++) {
         int idx  = IDX321   ( i, j,   k,   PS1+1, PS1+1 );
         int idxj = IDX321   ( i, j+1, k,   PS1+1, PS1+1 );
         int idxk = IDX321   ( i, j,   k+1, PS1+1, PS1+1 );
         int idxB = IDX321_BX( i, j,   k,   PS1,   PS1   );
         real Bx = ( Az[idxj] - Az[idx] - Ay[idxk] + Ay[idx] ) / dh;
         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[0][idxB] = Bx;
      }}}

//    Calculate By from vector potential
      for (int k=0; k<PS1;   k++) {
      for (int j=0; j<PS1+1; j++) {
      for (int i=0; i<PS1;   i++) {
         int idx  = IDX321   ( i,   j, k,   PS1+1, PS1+1 );
         int idxi = IDX321   ( i+1, j, k,   PS1+1, PS1+1 );
         int idxk = IDX321   ( i,   j, k+1, PS1+1, PS1+1 );
         int idxB = IDX321_BY( i,   j, k,   PS1,   PS1   );
         real By = ( Ax[idxk] - Ax[idx] - Az[idxi] + Az[idx] ) / dh;
         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[1][idxB] = By;
      }}}

//    Calculate Bz from vector potential
      for (int k=0; k<PS1+1; k++) {
      for (int j=0; j<PS1;   j++) {
      for (int i=0; i<PS1;   i++) {
         int idx  = IDX321   ( i,   j,   k, PS1+1, PS1+1 );
         int idxi = IDX321   ( i+1, j,   k, PS1+1, PS1+1 );
         int idxj = IDX321   ( i,   j+1, k, PS1+1, PS1+1 );
         int idxB = IDX321_BZ( i,   j,   k, PS1,   PS1   );
         real Bz = ( Ay[idxi] - Ay[idx] - Ax[idxj] + Ax[idx] ) / dh;
         amr->patch[ amr->MagSg[B_lv] ][B_lv][PID]->magnetic[2][idxB] = Bz;
      }}}

      delete [] Axf;
      delete [] Ayf;
      delete [] Azf;

   } // for (int PID=0; PID<amr->NPatchComma[B_lv][1]; PID++)

// Close the magnetic field file

#  ifdef SUPPORT_HDF5
   H5Fclose(mag_file_id);
#  endif

   delete [] Ax;
   delete [] Ay;
   delete [] Az;

   delete [] Axcoord;
   delete [] Aycoord;
   delete [] Azcoord;

// The patch data read once in VecPot_ReadPatches() (guarded by B_lv==0) is reused on
// every level, so it's only freed after the last level has been initialized.
   if ( B_lv == TOP_LEVEL ) {

      for (int p=0; p<NumBPatches; p++)  VecPot_FreePatch( BPatches[p] );

      delete [] BPatches;
      BPatches       = NULL;
      NumBPatches    = 0;
      MaxBPatchLevel = 0;

   }

   if ( MPI_Rank == 0 ) Aux_Message( stdout, "   Loading the magnetic field from the input file ... done\n" );

} // FUNCTION : MHD_Init_BField_ByVecPot_File



//-------------------------------------------------------------------------------------------------------
// Function    :  VecPot_Interp
// Description :  Use triangle-shaped cloud interpolation to interpolate
//                a vector potential from the input grid to the AMR grid
//
// Parameter   :  field  : Local patch of input vector potential, one component
//                xx     : coordinate along the x-axis
//                yy     : coordinate along the y-axis
//                zz     : coordinate along the z-axis
//                fdims  : size of the input vector potential patch
//                fbegin : index location of the input vector potential patch
//                grid   : coordinate arrays/spacing/lower-edge of the grid "field"
//                         lives on -- either the root coarse grid or a single
//                         refinement patch (see GridMeta)
//
// Return      :  vector potential component on AMR grid at (xx, yy, zz)
//-------------------------------------------------------------------------------------------------------
double VecPot_Interp( const double field[], const double xx, const double yy,
                                const double zz, const int fdims[], const int fbegin[],
                                const GridMeta &grid )
{

   // Indices into the coordinate vectors
   const int ii = (int)((xx-grid.xmin)/grid.dx);
   const int jj = (int)((yy-grid.ymin)/grid.dy);
   const int kk = (int)((zz-grid.zmin)/grid.dz);

// Indices into the local vector potential patch
   const int ib = ii - fbegin[0];
   const int jb = jj - fbegin[1];
   const int kb = kk - fbegin[2];

   double pot = 0.0;

   if ( ib == 0 || ib == fdims[0]-1 ||
        jb == 0 || jb == fdims[1]-1 ||
        kb == 0 || kb == fdims[2]-1 ) {
      Aux_Error( ERROR_INFO, "VecPot_Interp: An invalid index was entered!!\n" );
   }

   for (int i = -1; i <= 1; i++) { double dx = (xx-grid.xcoord[ii+i])/grid.dx;
   for (int j = -1; j <= 1; j++) { double dy = (yy-grid.ycoord[jj+j])/grid.dy;
   for (int k = -1; k <= 1; k++) { double dz = (zz-grid.zcoord[kk+k])/grid.dz;
      int idx = (ib+i)*fdims[2]*fdims[1] + (jb+j)*fdims[2] + (kb+k);
      pot += field[idx]*TSC_Weight(dx)*TSC_Weight(dy)*TSC_Weight(dz);
   }}}

   return pot;

} // FUNCTION : VecPot_Interp



//-------------------------------------------------------------------------------------------------------
// Function    :  VecPot_PointInPatch
// Description :  Determine whether (xx,yy,zz) is far enough inside the given refinement
//                patch's array to be safely interpolated by VecPot_Interp (which needs
//                a +/-1 cell margin around the nearest grid index). Unlike the coarse
//                grid (guaranteed by construction to cover the whole simulation domain
//                with margin, so an out-of-range index there is a fatal input error),
//                being outside this margin is the normal case for most query points --
//                refinement patches are small, local sub-boxes -- so this fails
//                gracefully (the caller falls back to the coarse-only value) instead of
//                raising an error.
//
// Parameter   :  patch : The refinement patch to test against
//                xx    : coordinate along the x-axis
//                yy    : coordinate along the y-axis
//                zz    : coordinate along the z-axis
//
// Return      :  true if (xx,yy,zz) can be safely interpolated on "patch"
//-------------------------------------------------------------------------------------------------------
bool VecPot_PointInPatch( const BFieldPatch &patch, const double xx, const double yy, const double zz )
{

   const int ii = (int)((xx-patch.xmin)/patch.dx);
   const int jj = (int)((yy-patch.ymin)/patch.dy);
   const int kk = (int)((zz-patch.zmin)/patch.dz);

   return ( ii >= 1 && ii <= patch.nx-2 &&
            jj >= 1 && jj <= patch.ny-2 &&
            kk >= 1 && kk <= patch.nz-2 );

} // FUNCTION : VecPot_PointInPatch

//-------------------------------------------------------------------------------------------------------
// Function    :  VecPot_EvalComponent
// Description :  Evaluate one component of the vector potential at (xx,yy,zz): the
//                finest refinement patch that contains the point, if any, or the
//                coarse (root-level) grid otherwise. When there are no patches at all
//                (NumBPatches == 0: a plain uniform-mesh file) this returns exactly
//                what VecPot_Interp on the coarse grid alone would, so the
//                uniform-grid-only case is unaffected.
//
//                Every patch's own data is already the final, fully-blended vector
//                potential (see BFieldPatch) -- not a correction on top of anything
//                coarser -- so there's no blending step here: this just searches from
//                the finest level (MaxBPatchLevel) down to 1 for the first (and, since
//                patches at the same level are assumed non-overlapping, only) patch
//                containing the point, and returns its value directly. Patches are
//                strictly nested, so the first match found searching from the finest
//                level down is always the most-refined patch covering the point.
//
// Parameter   :  coarse_field : Local buffer of the coarse grid's vector potential,
//                                one component (as read by VecPot_ReadField)
//                fdims        : size of the coarse_field buffer
//                fbegin       : index location of the coarse_field buffer
//                xx           : coordinate along the x-axis
//                yy           : coordinate along the y-axis
//                zz           : coordinate along the z-axis
//                comp         : which component: 0=x, 1=y, 2=z
//
// Return      :  vector potential component at (xx, yy, zz)
//-------------------------------------------------------------------------------------------------------
double VecPot_EvalComponent( const double coarse_field[], const int fdims[], const int fbegin[],
                             const double xx, const double yy, const double zz, const int comp )
{

   for (int lv=MaxBPatchLevel; lv>=1; lv--) {
   for (int p=0;  p<NumBPatches;        p++) {

      if ( BPatches[p].level == lv && VecPot_PointInPatch( BPatches[p], xx, yy, zz ) ) {

         const BFieldPatch &patch      = BPatches[p];
         const int          pfdims[3]  = { patch.nx, patch.ny, patch.nz };
         const int          pfbegin[3] = { 0, 0, 0 };
         const GridMeta     patch_grid = { patch.xcoord, patch.ycoord, patch.zcoord,
                                            patch.dx, patch.dy, patch.dz,
                                            patch.xmin, patch.ymin, patch.zmin };
         const double *comp_field = ( comp == 0 ) ? patch.Ax : ( comp == 1 ) ? patch.Ay : patch.Az;

         return VecPot_Interp( comp_field, xx, yy, zz, pfdims, pfbegin, patch_grid );

      }

   }}

   const GridMeta coarse_grid = { Axcoord, Aycoord, Azcoord, Adx, Ady, Adz, Axmin, Aymin, Azmin };

   return VecPot_Interp( coarse_field, xx, yy, zz, fdims, fbegin, coarse_grid );

} // FUNCTION : VecPot_EvalComponent

//-------------------------------------------------------------------------------------------------------
// Function    :  VecPot_FreePatch
// Description :  Free a single BFieldPatch entry's own arrays, as populated by
//                VecPot_ReadPatchGroup. Every patch -- top-level or nested -- is its own
//                independent entry in the flat BPatches array (see BFieldPatch), so,
//                unlike the old nested-group layout, this never needs to recurse.
//
// Parameter   :  patch : The patch to free
//
// Return      :  none
//-------------------------------------------------------------------------------------------------------
void VecPot_FreePatch( BFieldPatch &patch )
{

   delete [] patch.xcoord;
   delete [] patch.ycoord;
   delete [] patch.zcoord;
   delete [] patch.Ax;
   delete [] patch.Ay;
   delete [] patch.Az;

} // FUNCTION : VecPot_FreePatch

//-------------------------------------------------------------------------------------------------------
// Function    :  TSC_Weight
// Description :  Function to compute second-order interpolation on a uniform
//                grid to arbitrary points
//
// Parameter   :  x : Scaled coordinate for the TSC kernel
//
// Return      :  weight
//-------------------------------------------------------------------------------------------------------
double TSC_Weight( const double x )
{

   double weight;
   const double xx = fabs(x);

   if ( xx <= 0.5 ) {
      weight = 0.75 - SQR( xx );
   } else if ( xx >= 0.5 && xx <= 1.5 ) {
      weight = 0.5*SQR( 1.5-xx );
   } else {
      weight = 0.0;
   }

   return weight;

} // FUNCTION : TSC_Weight



#ifdef SUPPORT_HDF5
//-------------------------------------------------------------------------------------------------------
// Function    :  VecPot_ReadField 
// Description :  Read vector potential from a file 
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void VecPot_ReadField( hid_t mag_file_id, const int ibegin, const int jbegin,
                       const int kbegin, const int iend, const int jend,
                       const int kend, double Ax[], double Ay[], double Az[] )
{
   hid_t dataset, dataspace, memspace, dxfer_template;

   herr_t status;

   hsize_t start[3], stride[3], count[3], dims[3];

   int rank, ierr;

   rank = 3;

   start[0] = ibegin;
   start[1] = jbegin;
   start[2] = kbegin;

   stride[0] = 1;
   stride[1] = 1;
   stride[2] = 1;

   count[0] = iend-ibegin+1;
   count[1] = jend-jbegin+1;
   count[2] = kend-kbegin+1;

   dims[0] = count[0];
   dims[1] = count[1];
   dims[2] = count[2];

// Read Ax
   dataset = H5Dopen(mag_file_id, "magnetic_vector_potential_x", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                stride, count, NULL);
   memspace = H5Screate_simple(rank, dims, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, Ax);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load magnetic_vector_potential_x !!\n" );
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

// Read Ay
   dataset = H5Dopen(mag_file_id, "magnetic_vector_potential_y", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                 stride, count, NULL);
   memspace = H5Screate_simple(rank, dims, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                    H5P_DEFAULT, Ay);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load magnetic_vector_potential_y !!\n" );
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

// Read Az
   dataset = H5Dopen(mag_file_id, "magnetic_vector_potential_z", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                 stride, count, NULL);
   memspace = H5Screate_simple(rank, dims, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, Az);
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load magnetic_vector_potential_z !!\n" );
   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

   return;

} // FUNCTION : VecPot_ReadField

//-------------------------------------------------------------------------------------------------------
// Function    :  VecPot_ReadPatches
// Description :  Detect and read any locally-refined vector-potential patches (flat
//                HDF5 groups "patch_0000", "patch_0001", ..., optionally present in the
//                input file -- see BFieldPatch), populating the static BPatches array.
//
// Note        :  1. If the root-level "num_patches" attribute is absent -- a plain
//                   uniform-mesh file with no refinement patches at all -- NumBPatches
//                   is left at 0 and BPatches at NULL, so VecPot_EvalComponent falls
//                   back to exactly the coarse-grid-only behavior.
//                2. Each patch is small enough (by construction, a local sub-box of the
//                   full domain) that its full arrays -- not just a local hyperslab,
//                   unlike the coarse grid -- are read whole into every MPI rank.
//                3. Called once (guarded by B_lv==0 at the call site); the data is
//                   reused across all levels and freed only when B_lv==TOP_LEVEL.
//                4. Every patch -- however deeply nested -- is its own flat
//                   "patch_NNNN" group here, carrying no parent/child reference;
//                   VecPot_ReadPatchGroup reads each one's own data plus its "level"
//                   attribute (see BFieldPatch). MaxBPatchLevel (the highest level seen)
//                   is what lets VecPot_EvalComponent know how many levels to walk.
//
// Parameter   :  mag_file_id : Open HDF5 file identifier for the B-field input file
//
// Return      :  NumBPatches, MaxBPatchLevel, BPatches (static)
//-------------------------------------------------------------------------------------------------------
void VecPot_ReadPatches( hid_t mag_file_id )
{

   const htri_t has_patches = H5Aexists( mag_file_id, "num_patches" );

   if ( has_patches <= 0 ) {
      NumBPatches    = 0;
      MaxBPatchLevel = 0;
      BPatches       = NULL;
      return;
   }

   hid_t  attr_id = H5Aopen( mag_file_id, "num_patches", H5P_DEFAULT );
   herr_t status  = H5Aread( attr_id, H5T_NATIVE_INT, &NumBPatches );
   H5Aclose( attr_id );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load the \"num_patches\" attribute !!\n" );

   if ( NumBPatches <= 0 ) {
      NumBPatches    = 0;
      MaxBPatchLevel = 0;
      BPatches       = NULL;
      return;
   }

   BPatches = new BFieldPatch [ NumBPatches ];

   char group_name[32];

   MaxBPatchLevel = 0;

   for (int p=0; p<NumBPatches; p++) {

      snprintf( group_name, sizeof(group_name), "patch_%04d", p );

      hid_t group_id = H5Gopen( mag_file_id, group_name, H5P_DEFAULT );
      if ( group_id < 0 ) Aux_Error( ERROR_INFO, "Failed to open group \"%s\" !!\n", group_name );

      VecPot_ReadPatchGroup( group_id, BPatches[p] );

      H5Gclose( group_id );

      if ( BPatches[p].level > MaxBPatchLevel )  MaxBPatchLevel = BPatches[p].level;

   } // for (int p=0; p<NumBPatches; p++)

   return;

} // FUNCTION : VecPot_ReadPatches

//-------------------------------------------------------------------------------------------------------
// Function    :  VecPot_ReadPatchGroup
// Description :  Read a single flat "patch_NNNN" group's own data (coordinates,
//                Ax/Ay/Az, level) into "patch". Every patch -- however deeply nested --
//                is its own such group (see VecPot_ReadPatches), so this is never
//                called recursively; "patch.level" is this patch's own nesting depth
//                (1 for a patch nested directly inside the coarse grid, 2 for one
//                nested inside a level-1 patch, etc.) -- the only thing
//                VecPot_EvalComponent needs to pick out the finest patch covering a
//                given point.
//
// Parameter   :  group_id : Open HDF5 group identifier for this patch
//                patch    : Filled in place
//
// Return      :  "patch", filled in place
//-------------------------------------------------------------------------------------------------------
void VecPot_ReadPatchGroup( hid_t group_id, BFieldPatch &patch )
{

   hid_t   dataset, dataspace;
   hsize_t dims[3], maxdims[3];
   herr_t  status;

   char group_name[256];
   H5Iget_name( group_id, group_name, sizeof(group_name) );

// Dimensions, from the x/y/z coordinate datasets
   dataset   = H5Dopen( group_id, "x", H5P_DEFAULT );
   dataspace = H5Dget_space( dataset );
   H5Sget_simple_extent_dims( dataspace, dims, maxdims );
   patch.nx  = (int)dims[0];
   H5Sclose( dataspace );
   H5Dclose( dataset );

   dataset   = H5Dopen( group_id, "y", H5P_DEFAULT );
   dataspace = H5Dget_space( dataset );
   H5Sget_simple_extent_dims( dataspace, dims, maxdims );
   patch.ny  = (int)dims[0];
   H5Sclose( dataspace );
   H5Dclose( dataset );

   dataset   = H5Dopen( group_id, "z", H5P_DEFAULT );
   dataspace = H5Dget_space( dataset );
   H5Sget_simple_extent_dims( dataspace, dims, maxdims );
   patch.nz  = (int)dims[0];
   H5Sclose( dataspace );
   H5Dclose( dataset );

// Coordinate arrays
   patch.xcoord = new double [ patch.nx ];
   patch.ycoord = new double [ patch.ny ];
   patch.zcoord = new double [ patch.nz ];

   dataset = H5Dopen( group_id, "x", H5P_DEFAULT );
   status  = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, patch.xcoord );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/x\" !!\n", group_name );
   H5Dclose( dataset );

   dataset = H5Dopen( group_id, "y", H5P_DEFAULT );
   status  = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, patch.ycoord );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/y\" !!\n", group_name );
   H5Dclose( dataset );

   dataset = H5Dopen( group_id, "z", H5P_DEFAULT );
   status  = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, patch.zcoord );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/z\" !!\n", group_name );
   H5Dclose( dataset );

   patch.dx = patch.xcoord[1] - patch.xcoord[0];
   patch.dy = patch.ycoord[1] - patch.ycoord[0];
   patch.dz = patch.zcoord[1] - patch.zcoord[0];

   patch.xmin = patch.xcoord[0] - 0.5*patch.dx;
   patch.ymin = patch.ycoord[0] - 0.5*patch.dy;
   patch.zmin = patch.zcoord[0] - 0.5*patch.dz;

// Vector potential components, read whole (see Note 2 on VecPot_ReadPatches)
   const int npts = patch.nx*patch.ny*patch.nz;

   patch.Ax = new double [ npts ];
   patch.Ay = new double [ npts ];
   patch.Az = new double [ npts ];

   dataset = H5Dopen( group_id, "magnetic_vector_potential_x", H5P_DEFAULT );
   status  = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, patch.Ax );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/magnetic_vector_potential_x\" !!\n", group_name );
   H5Dclose( dataset );

   dataset = H5Dopen( group_id, "magnetic_vector_potential_y", H5P_DEFAULT );
   status  = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, patch.Ay );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/magnetic_vector_potential_y\" !!\n", group_name );
   H5Dclose( dataset );

   dataset = H5Dopen( group_id, "magnetic_vector_potential_z", H5P_DEFAULT );
   status  = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, patch.Az );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/magnetic_vector_potential_z\" !!\n", group_name );
   H5Dclose( dataset );

// level attribute -- this patch's own nesting depth (see BFieldPatch)
   hid_t level_attr = H5Aopen( group_id, "level", H5P_DEFAULT );
   status = H5Aread( level_attr, H5T_NATIVE_INT, &patch.level );
   H5Aclose( level_attr );
   if ( status < 0 ) Aux_Error( ERROR_INFO, "Failed to load \"%s/level\" !!\n", group_name );

} // FUNCTION : VecPot_ReadPatchGroup

#endif // #ifdef SUPPORT_HDF5
