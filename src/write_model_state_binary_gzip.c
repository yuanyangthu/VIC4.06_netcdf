#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
#include <zlib.h>

static char vcid[] = "$Id: write_model_state.c,v 4.2.2.4 2004/05/06 19:57:37 tbohn Exp $";

void write_model_state_binary_gzip(dist_prcp_struct    *prcp,
                                   global_param_struct *gp,
                                   int                  Nveg,
                                   int                  cellnum,
                                   outfiles_struct     *outfiles,
                                   soil_con_struct     *soil_con,
                                   char                 STILL_STORM,
                                   int                  DRY_TIME,
                                   dmy_struct           labeldate)
/*********************************************************************
  write_model_state      Keith Cherkauer           April 14, 2000

  This subroutine saves the model state at hour 0 of the date
  defined in the global control file using STATEDAY, STATEMONTH,
  and STATEYEAR.  The saved files can then be used to initialize
  the model to the same state as when the files were created.

  Soil moisture, soil thermal, and snowpack variables  are stored
  for each vegetation type and snow band.  However moisture variables
  from the distributed precipitation model are averaged so that the
  model is restarted with mu = 1.

  Modifications:
  04-10-03 Rewritten to handle updates to vicNl_def.h and to write
           the file as binary to minimize write time and differences
           with simulations started with the state file.         KAC
  04-10-03 Model is now restarted with the correct values for mu
           and LAST_STORM
  06-03-03 Modified to create ASCII as well as BINARY state file.  KAC
  09-05-2003 Modified to print space before dz_node for ASCII state
             file, this corrects a problem with state files created
             for models using the Cherkauer and Lettenmaier (1999) heat
             flux formulation.                                   KAC
  09-Oct-03 Added "\n" after mu in ASCII file, to jive with
            read_initial_model_state.                            TJB

*********************************************************************/
{
  extern option_struct options;
  extern global_param_struct global_param;

  double tmpval;
  double Nsum;
  int    veg;
  int    band;
  int    lidx;
  int    nidx;
  int    dist;
  int    Ndist;
  int    Nbands;
  int    byte, Nbytes;
  char   filename[MAXSTRING];

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;

  /* open state file */
  if (cellnum==-99999) {
    
    sprintf(filename,"%s_%04i%02i%02i.gz", gp->statename,
            labeldate.year, labeldate.month, labeldate.day);
    outfiles->gzstatefile = gzopen(filename,"wb");

    /* Write save state date information */
    gzwrite(outfiles->gzstatefile, (voidp *) &labeldate.year, sizeof(int));
    gzwrite(outfiles->gzstatefile, (voidp *) &labeldate.month, sizeof(int));
    gzwrite(outfiles->gzstatefile, (voidp *) &labeldate.day, sizeof(int));

    /* Write simulation flags */
    gzwrite(outfiles->gzstatefile, (voidp *) &options.Nlayer, sizeof(int));
    gzwrite(outfiles->gzstatefile, (voidp *) &options.Nnode, sizeof(int));
    
    return;
  }
    
  if(options.DIST_PRCP)
    Ndist = 2;
  else
    Ndist = 1;
  Nbands = options.SNOW_BAND;

  cell    = prcp->cell;
  veg_var = prcp->veg_var;
  snow    = prcp->snow;
  energy  = prcp->energy;

  /* write cell information */
  gzwrite(outfiles->gzstatefile, (voidp *) &cellnum, sizeof(int));
  gzwrite(outfiles->gzstatefile, (voidp *) &Nveg, sizeof(int));
  gzwrite(outfiles->gzstatefile, (voidp *) &Nbands, sizeof(int));
  // This stores the number of bytes from after this value to the end
  // of the line.  DO NOT CHANGE unless you have changed the values
  // written to the state file.
  // IF YOU EDIT THIS FILE: UPDATE THIS VALUE!
  Nbytes = ( options.Nnode * sizeof(double) // dz_node
             + sizeof(char) // STILL_STORM
             + sizeof(int) // DRY_TIME
             + Nveg * sizeof(double) // mu
             + Nveg * Nbands * 2 * sizeof(int) // veg & band
             + Nveg * Nbands * Ndist * options.Nlayer * sizeof(double) // soil moisture
             + Nveg * Nbands * Ndist * options.Nlayer * sizeof(double) // soil ice
             + (Nveg-1) * Nbands * sizeof(double) // dew
             + Nveg * Nbands * sizeof(int) // last_snow
             + Nveg * Nbands * sizeof(char) // MELTING
             + Nveg * Nbands * sizeof(double) * 9 // other snow parameters
             + Nveg * Nbands * options.Nnode * sizeof(double) ); // soil temperatures
  gzwrite(outfiles->gzstatefile, (voidp *) &Nbytes, sizeof(int));

  /* Write soil thermal node depths */
  Nsum = 0;
  for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
    gzwrite(outfiles->gzstatefile, (voidp *) &soil_con->dz_node[nidx], sizeof(double));
  }

  // Store distributed precipitation variables
  gzwrite(outfiles->gzstatefile, (voidp *) &STILL_STORM, sizeof(char));
  gzwrite(outfiles->gzstatefile, (voidp *) &DRY_TIME, sizeof(int));

  /* Output for all vegetation types */
  for ( veg = 0; veg <= Nveg; veg++ ) {

    // Store distributed precipitation fraction
    gzwrite(outfiles->gzstatefile, (voidp *) &prcp->mu[veg], sizeof(double));

    /* Output for all snow bands */
    for ( band = 0; band < Nbands; band++ ) {

      /* Write cell identification information */
      gzwrite(outfiles->gzstatefile, (voidp *) &veg, sizeof(int));
      gzwrite(outfiles->gzstatefile, (voidp *) &band, sizeof(int));

      for ( dist = 0; dist < Ndist; dist ++ ) {
        // Store both wet and dry fractions if using distributed precipitation

        /* Write total soil moisture */
        for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
          tmpval = cell[dist][veg][band].layer[lidx].moist;
          gzwrite(outfiles->gzstatefile, (voidp *) &tmpval, sizeof(double));
        }

        /* Write ice content */
        for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
          tmpval = cell[dist][veg][band].layer[lidx].ice;
          gzwrite(outfiles->gzstatefile, (voidp *) &tmpval, sizeof(double));
        }

        /* Write dew storage */
        if ( veg < Nveg ) {
          tmpval = veg_var[dist][veg][band].Wdew;
          gzwrite(outfiles->gzstatefile, (voidp *) &tmpval, sizeof(double));
        }
      }

      /* Write snow data */
      gzwrite(outfiles->gzstatefile, (voidp *) &snow[veg][band].last_snow, sizeof(int));
      gzwrite(outfiles->gzstatefile, (voidp *) &snow[veg][band].MELTING, sizeof(char));
      gzwrite(outfiles->gzstatefile, (voidp *) &snow[veg][band].coverage, sizeof(double));
      gzwrite(outfiles->gzstatefile, (voidp *) &snow[veg][band].swq, sizeof(double));
      gzwrite(outfiles->gzstatefile, (voidp *) &snow[veg][band].surf_temp, sizeof(double));
      gzwrite(outfiles->gzstatefile, (voidp *) &snow[veg][band].surf_water, sizeof(double));
      gzwrite(outfiles->gzstatefile, (voidp *) &snow[veg][band].pack_temp, sizeof(double));
      gzwrite(outfiles->gzstatefile, (voidp *) &snow[veg][band].pack_water, sizeof(double));
      gzwrite(outfiles->gzstatefile, (voidp *) &snow[veg][band].density, sizeof(double));
      gzwrite(outfiles->gzstatefile, (voidp *) &snow[veg][band].coldcontent, sizeof(double));
      gzwrite(outfiles->gzstatefile, (voidp *) &snow[veg][band].snow_canopy, sizeof(double));

      /* Write soil thermal node temperatures */
      for ( nidx = 0; nidx < options.Nnode; nidx++ )
        gzwrite(outfiles->gzstatefile, (voidp *) &energy[veg][band].T[nidx], sizeof(double));

    }
  }

  /* Force file to be written */
  // fflush(outfiles->statefile);

}
