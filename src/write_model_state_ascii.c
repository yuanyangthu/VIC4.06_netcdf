#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id: write_model_state.c,v 4.2.2.4 2004/05/06 19:57:37 tbohn Exp $";

void write_model_state(dist_prcp_struct    *prcp,
                       global_param_struct *gp,
                       int                  Nveg,
                       int                  cellnum,
                       outfiles_struct     *outfiles,
                       soil_con_struct     *soil_con,
                       char                 STILL_STORM,
                       int                  DRY_TIME,
                       dmy_struct           labeldate)
/* This is a wrapper to direct calls to the 4 functions */
{
  extern option_struct options;
  
  if (!options.BINARY_STATE_FILE) {
    if (!options.STATE_GZIP)
      write_model_state_ascii(prcp, gp, Nveg, cellnum, outfiles, soil_con, STILL_STORM, DRY_TIME, labeldate);
    else
      write_model_state_ascii_gzip(prcp, gp, Nveg, cellnum, outfiles, soil_con, STILL_STORM, DRY_TIME, labeldate);
  }
  else {
    if (!options.STATE_GZIP)
      write_model_state_binary(prcp, gp, Nveg, cellnum, outfiles, soil_con, STILL_STORM, DRY_TIME, labeldate);
    else
      write_model_state_binary_gzip(prcp, gp, Nveg, cellnum, outfiles, soil_con, STILL_STORM, DRY_TIME, labeldate);
  }
  
  return;  
}

void write_model_state_ascii(dist_prcp_struct    *prcp,
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
    
    sprintf(filename,"%s_%04d%02d%02d", gp->statename,
            labeldate.year, labeldate.month, labeldate.day);
    outfiles->statefile = open_file(filename,"w");

    fprintf(outfiles->statefile,"%d %d %d\n", labeldate.year, labeldate.month, labeldate.day);

    /* Write simulation flags */
    fprintf(outfiles->statefile,"%d %d\n", options.Nlayer, options.Nnode);
    
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
  fprintf( outfiles->statefile, "%i %i %i", cellnum, Nveg, Nbands );

  /* Write soil thermal node depths */
  Nsum = 0;
  for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
    fprintf( outfiles->statefile, " %f", soil_con->dz_node[nidx] );
  }
  fprintf( outfiles->statefile, "\n" );

  // Store distributed precipitation variables
  fprintf( outfiles->statefile, "%i %i\n", STILL_STORM, DRY_TIME );

  /* Output for all vegetation types */
  for ( veg = 0; veg <= Nveg; veg++ ) {

    // Store distributed precipitation fraction
    fprintf( outfiles->statefile, "%f\n", prcp->mu[veg] );

    /* Output for all snow bands */
    for ( band = 0; band < Nbands; band++ ) {

      /* Write cell identification information */
      fprintf( outfiles->statefile, "%i %i", veg, band );

      for ( dist = 0; dist < Ndist; dist ++ ) {
        // Store both wet and dry fractions if using distributed precipitation

        /* Write total soil moisture */
        for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
          tmpval = cell[dist][veg][band].layer[lidx].moist;
          fprintf( outfiles->statefile, " %f", tmpval );
        }

        /* Write ice content */
        for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
          tmpval = cell[dist][veg][band].layer[lidx].ice;
          fprintf( outfiles->statefile, " %f", tmpval );
        }

        /* Write dew storage */
        if ( veg < Nveg ) {
          tmpval = veg_var[dist][veg][band].Wdew;
          fprintf( outfiles->statefile, " %f", tmpval );
        }
      }

      /* Write snow data */
      fprintf( outfiles->statefile, " %i %i %f %f %f %f %f %f %f %f %f",
               snow[veg][band].last_snow, (int)snow[veg][band].MELTING,
               snow[veg][band].coverage, snow[veg][band].swq,
               snow[veg][band].surf_temp, snow[veg][band].surf_water,
               snow[veg][band].pack_temp, snow[veg][band].pack_water,
               snow[veg][band].density, snow[veg][band].coldcontent,
               snow[veg][band].snow_canopy );

      /* Write soil thermal node temperatures */
      for ( nidx = 0; nidx < options.Nnode; nidx++ )
        fprintf( outfiles->statefile, " %f", energy[veg][band].T[nidx] );

      fprintf( outfiles->statefile, "\n" );

    }
  }

  /* Force file to be written */
  //fflush(outfiles->statefile);

}
