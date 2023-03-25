#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id: write_model_state.c,v 4.2.2.4 2004/05/06 19:57:37 tbohn Exp $";

void write_model_state_binary(dist_prcp_struct    *prcp,
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
    
    sprintf(filename,"%s_%04i%02i%02i", gp->statename,
            labeldate.year, labeldate.month, labeldate.day);
    outfiles->statefile = open_file(filename,"wb");

    /* Write save state date information */
    fwrite( &labeldate.year, 1, sizeof(int), outfiles->statefile );
    fwrite( &labeldate.month, 1, sizeof(int), outfiles->statefile );
    fwrite( &labeldate.day, 1, sizeof(int), outfiles->statefile );

    /* Write simulation flags */
    fwrite( &options.Nlayer, 1, sizeof(int), outfiles->statefile );
    fwrite( &options.Nnode, 1, sizeof(int), outfiles->statefile );
    
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
  fwrite( &cellnum, 1, sizeof(int), outfiles->statefile );
  fwrite( &Nveg, 1, sizeof(int), outfiles->statefile );
  fwrite( &Nbands, 1, sizeof(int), outfiles->statefile );
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
  fwrite( &Nbytes, 1, sizeof(int), outfiles->statefile );

  /* Write soil thermal node depths */
  Nsum = 0;
  for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
    fwrite( &soil_con->dz_node[nidx], 1, sizeof(double),
            outfiles->statefile );
  }

  // Store distributed precipitation variables
  fwrite( &STILL_STORM, 1, sizeof(char), outfiles->statefile );
  fwrite( &DRY_TIME, 1, sizeof(int), outfiles->statefile );

  /* Output for all vegetation types */
  for ( veg = 0; veg <= Nveg; veg++ ) {

    // Store distributed precipitation fraction
    fwrite( &prcp->mu[veg], 1, sizeof(double), outfiles->statefile );

    /* Output for all snow bands */
    for ( band = 0; band < Nbands; band++ ) {

      /* Write cell identification information */
      fwrite( &veg, 1, sizeof(int), outfiles->statefile );
      fwrite( &band, 1, sizeof(int), outfiles->statefile );

      for ( dist = 0; dist < Ndist; dist ++ ) {
        // Store both wet and dry fractions if using distributed precipitation

        /* Write total soil moisture */
        for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
          tmpval = cell[dist][veg][band].layer[lidx].moist;
          fwrite( &tmpval, 1, sizeof(double), outfiles->statefile );
        }

        /* Write ice content */
        for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
          tmpval = cell[dist][veg][band].layer[lidx].ice;
          fwrite( &tmpval, 1, sizeof(double), outfiles->statefile );
        }

        /* Write dew storage */
        if ( veg < Nveg ) {
          tmpval = veg_var[dist][veg][band].Wdew;
          fwrite( &tmpval, 1, sizeof(double), outfiles->statefile );
        }
      }

      /* Write snow data */
      fwrite( &snow[veg][band].last_snow, 1, sizeof(int), outfiles->statefile );
      fwrite( &snow[veg][band].MELTING, 1, sizeof(char), outfiles->statefile );
      fwrite( &snow[veg][band].coverage, 1, sizeof(double), outfiles->statefile );
      fwrite( &snow[veg][band].swq, 1, sizeof(double), outfiles->statefile );
      fwrite( &snow[veg][band].surf_temp, 1, sizeof(double), outfiles->statefile );
      fwrite( &snow[veg][band].surf_water, 1, sizeof(double), outfiles->statefile );
      fwrite( &snow[veg][band].pack_temp, 1, sizeof(double), outfiles->statefile );
      fwrite( &snow[veg][band].pack_water, 1, sizeof(double), outfiles->statefile );
      fwrite( &snow[veg][band].density, 1, sizeof(double), outfiles->statefile );
      fwrite( &snow[veg][band].coldcontent, 1, sizeof(double), outfiles->statefile );
      fwrite( &snow[veg][band].snow_canopy, 1, sizeof(double), outfiles->statefile );

      /* Write soil thermal node temperatures */
      for ( nidx = 0; nidx < options.Nnode; nidx++ )
        fwrite( &energy[veg][band].T[nidx], 1, sizeof(double),
                outfiles->statefile );

    }
  }

  /* Force file to be written */
  //fflush(outfiles->statefile);

}
