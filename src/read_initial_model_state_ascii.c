#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id: read_initial_model_state.c,v 4.3.2.4 2004/05/11 19:39:56 tbohn Exp $";

void read_initial_model_state_ascii(FILE                *statefile,
                                    dist_prcp_struct    *prcp,
                                    global_param_struct *gp,
                                    int                  Nveg,
                                    int                  Nbands,
                                    int                  cellnum,
                                    soil_con_struct     *soil_con,
                                    int                  Ndist,
                                    char                *init_STILL_STORM,
                                    int                 *init_DRY_TIME)
/*********************************************************************
  read_initial_model_state   Keith Cherkauer         April 14, 2000

  This subroutine initializes the model state at hour 0 of the date
  defined in the given state file.

  Soil moisture, soil thermal, and snowpack variables  are stored
  for each vegetation type and snow band.  However moisture variables
  from the distributed precipitation model are averaged so that the
  model is restarted with mu = 1.

  Modifications:
  04-10-03 Rewritten to handle updates to vicNl_def.h and to write
           the file as binary to minimize write time and differences
           with simulations started with the state file.         KAC
  04-10-03 Modified to read storm parameters from the state file.  KAC
  06-03-03 Modified to read ASCII as well as BINARY state file.  KAC
  03-Oct-03 Modified to loop over tmp_Nveg and tmp_Nband when searching
            for desired cellnum in ASCII file, rather than over Nveg
            and Nbands.  As we skip over other records in the state
            file while searching for the desired record, the loop
            must parse each undesired record differently, according
            to how many veg classes and snow bands exist in the
            record (tmp_Nveg and tmp_Nband, respectively), rather
            than the number of veg classes and snow bands in the
            desired record (Nveg and Nbands, respectively).       TJB
  09-Oct-03 Modified to skip an extra line (for init_STILL_STORM
            and init_DRY_TIME) when searching for desired cellnum in
            ASCII file.                                           TJB
  11-May-04 (Port from 4.1.0) Added check to verify that the sum of
            the defined nodes equals the damping depth.                TJB

*********************************************************************/
{
  extern option_struct options;

  char   tmpstr[MAXSTRING];
  char   ErrStr[MAXSTRING];
  char   tmpchar;
  double tmpval;
  double Nsum;
  double depth_node[MAX_NODES];
  double sum;
  int    veg, iveg;
  int    band, iband;
  int    lidx;
  int    nidx;
  int    dist;
  int    tmp_cellnum;
  int    tmp_Nveg;
  int    tmp_Nband;
  int    tmp_char;
  int    byte, Nbytes;

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;

  cell    = prcp->cell;
  veg_var = prcp->veg_var;
  snow    = prcp->snow;
  energy  = prcp->energy;

  /* read cell information */
  //fscanf(statefile, "%i %i %i", &tmp_cellnum, &tmp_Nveg, &tmp_Nband);
  fscanf( statefile, "%i %i %i", &tmp_cellnum, &tmp_Nveg, &tmp_Nband );

  // Skip over unused cell information
  while ( tmp_cellnum != cellnum && !feof(statefile) ) {
      // skip rest of current cells info
      fgets(tmpstr, MAXSTRING, statefile); // skip remainder of this line (soil thermal node depths)
      fgets(tmpstr, MAXSTRING, statefile); // skip init_STILL_STORM, init_DRY_TIME
      for ( veg = 0; veg <= tmp_Nveg; veg++ ) {
        fgets(tmpstr, MAXSTRING, statefile); // skip dist precip info
        for ( band = 0; band < tmp_Nband; band++ )
          fgets(tmpstr, MAXSTRING, statefile); // skip veg/snowband info
      }
      // read info for next cell
      fscanf( statefile, "%i %i %i", &tmp_cellnum, &tmp_Nveg, &tmp_Nband );
  }

  if ( feof(statefile) ) {
    sprintf(ErrStr, "Requested grid cell (%i) is not in the model state file.",
            cellnum);
    nrerror(ErrStr);
  }

  if ( tmp_Nveg != Nveg ) {
    sprintf(ErrStr,"The number of vegetation types in cell %i (%i) does not equal that defined in vegetation parameter file (%i).  Check your input files.", cellnum, tmp_Nveg, Nveg);
    nrerror(ErrStr);
  }
  
  if ( tmp_Nband != Nbands ) {
    sprintf(ErrStr,"The number of snow bands in cell %i (%i) does not equal that defined in the snow band file (%i).  Check your input files.", cellnum, tmp_Nband, Nbands);
    nrerror(ErrStr);
  }
  /*
  if ( tmp_Nband < soil_con->n_bands ) {
    sprintf(ErrStr,"The number of snow bands in cell %i (%i) is less than the number of non-zero bands defined in the snow band file (%i).  Check your input files.",
            cellnum, tmp_Nband, soil_con->n_bands);
    nrerror(ErrStr);
  }
  */

  /* Read soil thermal node depths */
  sum = 0;
  for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
    fscanf( statefile, "%lf", &soil_con->dz_node[nidx] );
    sum += soil_con->dz_node[nidx];
  }
  if ( options.Nnode == 1 ) soil_con->dz_node[0] = 0;
  sum -= 0.5 * soil_con->dz_node[0];
  sum -= 0.5 * soil_con->dz_node[options.Nnode-1];
  if ( abs( sum - soil_con->dp ) > SMALL ) {
    fprintf( stderr, "WARNING: Sum of soil nodes (%f) exceeds defined damping depth (%f).  Resetting damping depth.\n", sum,
soil_con->dp );
    soil_con->dp = sum;
  }

  // read storm parameters
  fscanf( statefile, "%i %i", &tmp_char, init_DRY_TIME );
  (*init_STILL_STORM) = (char)tmp_char;

  /* Input for all vegetation types */
  for ( veg = 0; veg <= Nveg; veg++ ) {

    // read distributed precipitation fraction
    fscanf( statefile, "%lf", &prcp->mu[veg] );

    /* Input for all snow bands */
    for ( band = 0; band < tmp_Nband; band++ ) {

      /* skip unused bands to make it compatible with old state files */
      if (band>=Nbands) {
        if (band==Nbands) fgets(tmpstr, MAXSTRING, statefile);
        fgets(tmpstr, MAXSTRING, statefile);
        //fprintf(stderr, "skipping band %d\n", band);
        continue;
      }
      /* Read cell identification information */
      if ( fscanf(statefile,"%i %i", &iveg, &iband) == EOF )
        nrerror("End of model state file found unexpectedly");
      if ( iveg != veg || iband != band ) {
        fprintf(stderr,"The vegetation and snow band indices in the model state file (veg = %i, band = %i) do not match those currently requested (veg = %i , band = %i).  Model state file must be stored with variables for all vegetation indexed by variables for all snow bands.", iveg, iband, veg, band);
        nrerror(ErrStr);
      }

      // Read both wet and dry fractions if using distributed precipitation
      for ( dist = 0; dist < Ndist; dist ++ ) {

        /* Read total soil moisture */
        for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
          if ( fscanf(statefile," %lf",
                      &cell[dist][veg][band].layer[lidx].moist) == EOF )
            nrerror("End of model state file found unexpectedly");
        }

        /* Read ice content */
        for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
          if ( fscanf(statefile," %lf",
                      &cell[dist][veg][band].layer[lidx].ice) == EOF )
            nrerror("End of model state file found unexpectedly");
        }

        /* Read dew storage */
        if ( veg < Nveg ) {
          if ( fscanf(statefile," %lf", &veg_var[dist][veg][band].Wdew) == EOF )
            nrerror("End of model state file found unexpectedly");
        }
      }

      /* Read snow data */
      if ( fscanf(statefile," %i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                  &snow[veg][band].last_snow, &tmp_char,
                  &snow[veg][band].coverage, &snow[veg][band].swq,
                  &snow[veg][band].surf_temp, &snow[veg][band].surf_water,
                  &snow[veg][band].pack_temp, &snow[veg][band].pack_water,
                  &snow[veg][band].density, &snow[veg][band].coldcontent,
                  &snow[veg][band].snow_canopy)
           == EOF )
        nrerror("End of model state file found unexpectedly");
      snow[veg][band].MELTING = (char)tmp_char;
      if(snow[veg][band].density > 0.)
        snow[veg][band].depth = 1000. * snow[veg][band].swq
          / snow[veg][band].density;

      /* Read soil thermal node temperatures */
      for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
        if ( fscanf(statefile," %lf", &energy[veg][band].T[nidx]) == EOF )
          nrerror("End of model state file found unexpectedly");

      }
    }
  }
}
