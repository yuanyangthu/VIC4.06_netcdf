#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
#include <global.h>

static char vcid[] = "$Id: vicNl.c,v 4.2.2.4 2004/06/17 21:47:47 tbohn Exp $";

/** Main Program **/

int main(int argc, char *argv[])
/**********************************************************************
        vicNl.c                Dag Lohmann                January 1996

  This program controls file I/O and variable initialization as well as
  being the primary driver for the model.

  For details about variables, input files and subroutines check:
        http://ce.washington.edu/~hydro/Lettenmaier/Models/VIC/VIC_home.html

  UNITS: unless otherwise marked:
         all water balance components are in mm
         all energy balance components are in mks
         depths, and lengths are in m

  modifications:
  1997-98 Model was updated from simple 2 layer water balance to 
          an extension of the full energy and water balance 3 layer
          model.                                                  KAC
  03-12-03 Modifed to add AboveTreeLine to soil_con_struct so that
           the model can make use of the computed treeline.     KAC
  04-10-03 Modified to initialize storm parameters using the state
           file.                                                KAC
  04-10-03 Modified to start the model by skipping records until the
           state file date is found.  This replaces the previous method
           of modifying the global file start date, which can change 
           the interpolation of atmospheric forcing data.        KAC
  04-15-03 Modified to store wet and dry fractions when intializing 
           water balance storage.  This accounts for changes in model
           state initialization, which now stores wet and dry fractions
           rather than just averagedvalues.                      KAC
  29-Oct-03 Modified the version display banner to print the version
            string defined in global.h.                                TJB
  16-Jun-04 Modified to pass soil_con.avgJulyAirTemp to
            initialize_atmos().                                        TJB

**********************************************************************/
{

  extern veg_lib_struct *veg_lib;
  extern option_struct options;
  extern Error_struct Error;
  extern global_param_struct global_param;
  extern param_set_struct param_set;

  /** Variable Declarations **/

  char                     LASTREC;
  char                     MODEL_DONE;
  char                     init_STILL_STORM, *init_STILL_STORM_all;
  int                      rec, i, j;
  int                      veg;
  int                      dist;
  int                      band;
  int                      Ndist;
  int                      Nveg_type;
  int                      index;
  int                      init_DRY_TIME, *init_DRY_TIME_all;
  int                      RUN_MODEL;
  int                      cell_cnt, dummy;
  int                      startrec;
  double                   storage;
  double                   veg_fract;
  double                   band_fract;
  dmy_struct              *dmy;
  atmos_data_struct       *atmos;
  veg_con_struct          *veg_con, *veg_con2, **veg_con_all;
  soil_con_struct          soil_con, *soil_con_all;
  dist_prcp_struct         prcp, *prcp_all; /* stores information about distributed 
                                    precipitation */
  veg_lib_struct          *veg_lib_tmp;
  FILE                    *fp;
  
  filenames_struct         filenames;
  filenames_struct         builtnames;
  infiles_struct           infiles;
  outfiles_struct          outfiles;
  
  float                  **nldas_forcing = (float **)NULL; /* stores all the forcing data */
  char                     linebuffer[MAXSTRING];
  
  dist_prcp_struct **prcp_states; /* to store monthly model states */
  dmy_struct        *dmy_states = NULL;  /* dates to dump model states */
  int                i_state, n_states=0;    /* number of dates for state dumping */
  char             **init_STILL_STORM_states;
  int              **init_DRY_TIME_states;
  storm_track_struct tmp_storm;
  
  tmp_storm.still_storm = FALSE;
  tmp_storm.dry_time    = -999;
  
  /** Read Model Options **/
  initialize_global();
  filenames = cmd_proc(argc, argv);

  /** Read Global Control File **/
  infiles.globalparam = open_file(filenames.global,"r");
  global_param = get_global_param(&filenames, infiles.globalparam);

  /** Check and Open Files **/
  check_files(&infiles, &filenames);

  /** Read Vegetation Library File **/
  veg_lib = read_veglib(infiles.veglib,&Nveg_type);

  /** Initialize Parameters **/
  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;

  /** Make Date Data Structure **/
  dmy      = make_dmy(&global_param);

  if ( options.INIT_STATE ) {
    if (!options.STATE_GZIP)
      infiles.statefile = check_state_file(filenames.init_state, dmy,
                                           &global_param, options.Nlayer,
                                           options.Nnode, &startrec);
    else
      infiles.gzstatefile = check_state_file_gzip(filenames.init_state, dmy,
                                                &global_param, options.Nlayer,
                                                options.Nnode, &startrec);
  }
  else
    startrec = 0;
  
  //fprintf(stderr, "Start reading input forcing ...\n");
  /** Gridded Iutput, read the whole forcing file in once **/
  if (options.GRID_INPUT) {
      if (param_set.FORCE_FORMAT[0] == BINARY)
          nldas_forcing = read_once_nldas_forcing(&filenames);
      else
          nldas_forcing = read_once_netcdf_forcing(&filenames);
  }

  fprintf(stderr, "Initializing output ...\n");
  /** Gridded Output: initialize **/
  if (options.GRID_OUTPUT) {
      if (options.BINARY_OUTPUT)
          write_data_nldas_grid((out_data_struct *)NULL, &filenames, &global_param, -9999, global_param.nrecs, -1);
      else
          write_data_netcdf((out_data_struct *)NULL, &filenames, &global_param, -9999, global_param.nrecs, -1);
  }
  
  /** find out how many monthly dates for state dumping **/
  if ( strcmp( global_param.statename, "NONE" ) != 0 && strcmp( global_param.statefreq, "MONTH" ) == 0) {
      
      dmy_states = (dmy_struct *) calloc(MAX_STATES, sizeof(dmy_struct));
      n_states = 0;
      for (rec=0; rec<global_param.nrecs-1; rec++) {  /* skip the first and last recs (starting/ending dates) */
          if (dmy[rec+1].day==1 && n_states<MAX_STATES) {
              dmy_states[n_states] = dmy[rec+1];
              n_states++;
          }
      }
      fprintf(stderr, "%d of monthly model states + 1 state at the end of the simulation to save/dump ...\n", n_states);
      
      prcp_states = (dist_prcp_struct **) calloc(n_states, sizeof(dist_prcp_struct *));
      init_STILL_STORM_states = (char **) calloc(n_states, sizeof(char *));
      init_DRY_TIME_states    = (int **)  calloc(n_states, sizeof(int *));
      
      for (i_state=0; i_state<n_states; i_state++) {
          prcp_states[i_state] = (dist_prcp_struct *) calloc(global_param.ncells, sizeof(dist_prcp_struct));
          init_STILL_STORM_states[i_state] = (char *) calloc(global_param.ncells, sizeof(char));
          init_DRY_TIME_states[i_state]    = (int *)  calloc(global_param.ncells, sizeof(int));
      }
  }
  
  /** read in all soil params for all cells at onece **/
  
  soil_con_all = (soil_con_struct *) calloc(global_param.ncells, sizeof(soil_con_struct));
  veg_con_all  = (veg_con_struct **) calloc(global_param.ncells, sizeof(veg_con_struct *));
  prcp_all     = (dist_prcp_struct *) calloc(global_param.ncells, sizeof(dist_prcp_struct));
  init_STILL_STORM_all = (char *) calloc(global_param.ncells, sizeof(char));
  init_DRY_TIME_all    = (int *) calloc(global_param.ncells, sizeof(int));
  
  fprintf(stderr, "Reading soil/veg/snowband/init ...\n");
  cell_cnt=0;
  while ((fscanf(infiles.soilparam, "%d", &flag))!=EOF) {
      if (flag) {
          /** read soil parameters **/
          soil_con_all[cell_cnt] = read_soilparam(infiles.soilparam, TRUE);
          /** read vegetation parameters **/
          veg_con_all[cell_cnt]  = read_vegparam(infiles.vegparam, soil_con_all[cell_cnt].gridcel, Nveg_type);
          calc_root_fractions(veg_con_all[cell_cnt], &soil_con_all[cell_cnt]);
          /** read snow band parameters **/
          read_snowband(infiles.snowband,soil_con_all[cell_cnt].gridcel,
                    (double)soil_con_all[cell_cnt].elevation, &soil_con_all[cell_cnt].Tfactor, 
                    &soil_con_all[cell_cnt].Pfactor, &soil_con_all[cell_cnt].AreaFract, 
                    &soil_con_all[cell_cnt].AboveTreeLine);
          /** Make Precipitation Distribution Control Structure **/
          prcp_all[cell_cnt] = make_dist_prcp(veg_con_all[cell_cnt][0].vegetat_type_num, &options.Nnode);
          if ( strcmp( global_param.statename, "NONE" ) != 0 && strcmp( global_param.statefreq, "MONTH" ) == 0) {
              for (i_state=0; i_state<n_states; i_state++)
                  prcp_states[i_state][cell_cnt] = make_dist_prcp(veg_con_all[cell_cnt][0].vegetat_type_num, &options.Nnode);
          }
          
          /** initialize model state **/
          initialize_model_state(&prcp_all[cell_cnt], dmy[0], soil_con_all[cell_cnt].avg_temp, 
                             &global_param, infiles, soil_con_all[cell_cnt].gridcel, 
                             veg_con_all[cell_cnt][0].vegetat_type_num, options.Nnode, 
                             Ndist, &soil_con_all[cell_cnt], veg_con_all[cell_cnt], &init_STILL_STORM_all[cell_cnt],
                             &init_DRY_TIME_all[cell_cnt]);
          cell_cnt++;
      }
      else fgets(linebuffer, MAXSTRING, infiles.soilparam);
  }
  
  if (options.TRUE_LATLON) {
      fprintf(stderr, "Reading true lat/lon ...\n");
      fp = open_file(filenames.true_latlon, "r");
      fgets(linebuffer, MAXSTRING, fp);
      cell_cnt = 0;
      while (!feof(fp)) {
          sscanf(linebuffer, "%d %d %f %f", &flag, &dummy, &soil_con_all[cell_cnt].lat, &soil_con_all[cell_cnt].lng);
          soil_con_all[cell_cnt].time_zone_lng = round(soil_con_all[cell_cnt].lng/15.0) * 15.0;
          if (flag) cell_cnt++;
          fgets(linebuffer, MAXSTRING, fp);
      }
      fclose(fp);
  }
  
  /** Build Gridded Filenames, and Open **/
  builtnames = make_in_and_outfiles(&infiles, &filenames, &soil_con_all[0],
               &outfiles, &global_param);
  /** Update Error Handling Structure **/
  Error.outfp = outfiles;
  Error.infp = infiles;

  fprintf(stderr, "Start parallel cell loop ...\n");
  /************************************
    Run Model for all Active Grid Cells
    ************************************/
  #pragma omp parallel for private(soil_con, veg_con, veg_lib_tmp, j, prcp, init_STILL_STORM, init_DRY_TIME, atmos, rec, LASTREC, storage, veg_fract, band_fract, veg, band, index, dist, i_state, tmp_storm)
  for (cell_cnt=0; cell_cnt<global_param.ncells; cell_cnt++) {
    
      soil_con = soil_con_all[cell_cnt];
      veg_con  = veg_con_all[cell_cnt];
      
      veg_lib_tmp = (veg_lib_struct *)calloc(Nveg_type,sizeof(veg_lib_struct));
      for ( veg = 0; veg < Nveg_type; veg++ ) {
          veg_lib_tmp[veg].overstory   = veg_lib[veg].overstory;
          veg_lib_tmp[veg].rad_atten   = veg_lib[veg].rad_atten;
          veg_lib_tmp[veg].rarc        = veg_lib[veg].rarc;
          veg_lib_tmp[veg].rmin        = veg_lib[veg].rmin;
          veg_lib_tmp[veg].trunk_ratio = veg_lib[veg].trunk_ratio;
          veg_lib_tmp[veg].wind_atten  = veg_lib[veg].wind_atten;
          veg_lib_tmp[veg].wind_h      = veg_lib[veg].wind_h;
          veg_lib_tmp[veg].RGL         = veg_lib[veg].RGL;
          veg_lib_tmp[veg].veg_class   = veg_lib[veg].veg_class;
          for (j=0; j<12; j++) {
              veg_lib_tmp[veg].LAI[j]          = veg_lib[veg].LAI[j];
              veg_lib_tmp[veg].Wdmax[j]        = veg_lib[veg].Wdmax[j];
              veg_lib_tmp[veg].albedo[j]       = veg_lib[veg].albedo[j];
              veg_lib_tmp[veg].displacement[j] = veg_lib[veg].displacement[j];
              veg_lib_tmp[veg].emissivity[j]   = veg_lib[veg].emissivity[j];
              veg_lib_tmp[veg].roughness[j]    = veg_lib[veg].roughness[j];
          }
      }
      
      /** when GLOBAL_LAI is turned, veg_lib needs to be updated per cell **/
      if ( options.GLOBAL_LAI )
          for ( veg = 0; veg < veg_con[0].vegetat_type_num; veg++ )
              for (j=0; j<12; j++) {
                  veg_lib_tmp[veg_con[veg].veg_class].LAI[j]   = veg_con[veg].LAI[j];
                  veg_lib_tmp[veg_con[veg].veg_class].Wdmax[j] = veg_con[veg].Wdmax[j];
              }
      
      prcp     = prcp_all[cell_cnt];
      init_STILL_STORM = init_STILL_STORM_all[cell_cnt];
      init_DRY_TIME    = init_DRY_TIME_all[cell_cnt];
      
      /** allocate memory for the atmos_data_struct **/
      alloc_atmos(global_param.nrecs, &atmos);
      /**************************************************
         Initialize Meteological Forcing Values That
         Have not Been Specifically Set
       **************************************************/

      initialize_atmos(atmos, dmy, infiles.forcing, 
                       (double)soil_con.time_zone_lng, (double)soil_con.lng,
                       (double)soil_con.lat, soil_con.elevation,
                       soil_con.annual_prec, global_param.wind_h, 
                       soil_con.rough, soil_con.avgJulyAirTemp, 
                       soil_con.Tfactor,
                       nldas_forcing, cell_cnt, /* Gridded Input */
                       soil_con.AboveTreeLine); 



      /***************************************************
        Intialize Moisture and Energy Balance Error Checks
        --- As of 4/15/03 this does not properly initialize
            storage from bands above treeline, when the model 
            state is restored from a file.  This can lead to 
            water balance errors in the initial time step but 
            does not impact the actual simulation.  It will
            be addressed in the next release version.  KAC
        ***************************************************/
      storage = 0.;
      for ( veg = 0; veg <= veg_con[0].vegetat_type_num; veg++ ) {
        if ( veg < veg_con[0].vegetat_type_num ) veg_fract = veg_con[veg].Cv;
        else veg_fract = ( 1.0 - veg_con[0].Cv_sum );
        for ( band = 0; band < options.SNOW_BAND; band++ ) {
          band_fract = soil_con.AreaFract[band];
          if ( veg_fract > SMALL && band_fract > SMALL ) {
            for(index=0;index<options.Nlayer;index++)
              for ( dist = 0; dist < Ndist; dist ++ )
                storage += prcp.cell[dist][veg][band].layer[index].moist 
                  * veg_fract * band_fract;
            storage += prcp.snow[veg][band].swq * 1000. * veg_fract 
              * band_fract;
            if ( veg != veg_con[0].vegetat_type_num ) {
              for ( dist = 0; dist < Ndist; dist ++ ) 
                storage += prcp.veg_var[dist][veg][band].Wdew 
                  * veg_fract * band_fract;
              storage += prcp.snow[veg][band].snow_canopy * 1000. 
                * veg_fract * band_fract;
            }
          }
        }
      }
      calc_water_balance_error(-global_param.nrecs,0.,0.,storage);
      calc_energy_balance_error(-global_param.nrecs,0.,0.,0.,0.,0.);

      /******************************************
        Run Model in Grid Cell for all Time Steps
        ******************************************/

      i_state = 0;
      for ( rec = startrec ; rec < global_param.nrecs; rec++ ) {

        if ( rec == global_param.nrecs - 1 ) LASTREC = TRUE;
        else LASTREC = FALSE;

        tmp_storm = dist_prec( &atmos[rec], &prcp, &soil_con, veg_con,
                               dmy, &global_param, &outfiles, rec,
                               LASTREC, init_STILL_STORM, init_DRY_TIME, cell_cnt, veg_lib_tmp);
        init_DRY_TIME = -999;
        
        if (strcmp( global_param.statename, "NONE" )!=0 && strcmp( global_param.statefreq, "MONTH" )==0
            && rec<global_param.nrecs-1 && dmy[rec+1].day==1 && i_state<MAX_STATES) {
            copy_dist_prcp(&prcp, &prcp_states[i_state][cell_cnt], veg_con[0].vegetat_type_num);
            init_STILL_STORM_states[i_state][cell_cnt] = tmp_storm.still_storm;
            init_DRY_TIME_states[i_state][cell_cnt]    = tmp_storm.dry_time;
            i_state++;
        }

      }        /* End Rec Loop */

      init_STILL_STORM_all[cell_cnt] = tmp_storm.still_storm;
      init_DRY_TIME_all[cell_cnt]    = tmp_storm.dry_time;
      
      free(veg_lib_tmp);
      free_atmos(global_param.nrecs, &atmos);

  }         /* End Grid Loop */
  
  close_files(&infiles,&outfiles,&builtnames); 

  /** Gridded Output: write to file **/
  if (options.GRID_OUTPUT) {
      if (options.BINARY_OUTPUT)
          write_data_nldas_grid((out_data_struct *)NULL, &filenames, &global_param, -9999, global_param.nrecs, -2);
      else
          write_data_netcdf((out_data_struct *)NULL, &filenames, &global_param, -9999, global_param.nrecs, -2);
  }

  /** write model state file **/
  if ( strcmp( global_param.statename, "NONE" ) != 0 ) {
      
      fprintf(stderr, "Writing statefile ...\n");
      write_model_state(&prcp_all[0], &global_param, veg_con_all[0][0].vegetat_type_num,
                        -99999, &outfiles, &soil_con_all[0],
                        init_STILL_STORM_all[0], init_DRY_TIME_all[0], dmy[global_param.nrecs]);
      
      for (cell_cnt=0; cell_cnt<global_param.ncells; cell_cnt++)
          write_model_state(&prcp_all[cell_cnt], &global_param, veg_con_all[cell_cnt][0].vegetat_type_num,
                            soil_con_all[cell_cnt].gridcel, &outfiles, &soil_con_all[cell_cnt],
                            init_STILL_STORM_all[cell_cnt], init_DRY_TIME_all[cell_cnt], dmy[global_param.nrecs]);
      
      if (!options.STATE_GZIP) fclose(outfiles.statefile);
      else                     gzclose(outfiles.gzstatefile);
      
      /** write monthly states **/
      if (strcmp( global_param.statefreq, "MONTH" )==0) {
          for (i_state=0; i_state<n_states; i_state++) {
              
              write_model_state(&prcp_states[i_state][0], &global_param, veg_con_all[0][0].vegetat_type_num,
                                -99999, &outfiles, &soil_con_all[0],
                                init_STILL_STORM_states[i_state][0], init_DRY_TIME_states[i_state][0], dmy_states[i_state]);
              for (cell_cnt=0; cell_cnt<global_param.ncells; cell_cnt++)
                  write_model_state(&prcp_states[i_state][cell_cnt], &global_param, veg_con_all[cell_cnt][0].vegetat_type_num,
                            soil_con_all[cell_cnt].gridcel, &outfiles, &soil_con_all[cell_cnt],
                            init_STILL_STORM_states[i_state][cell_cnt], init_DRY_TIME_states[i_state][cell_cnt], dmy_states[i_state]);
      
              if (!options.STATE_GZIP) fclose(outfiles.statefile);
              else                     gzclose(outfiles.gzstatefile);
              
          }
      }
      
  }
      
  /** cleanup **/  
  if (strcmp( global_param.statename, "NONE" )!=0 && strcmp( global_param.statefreq, "MONTH" )==0) {
      for (i_state=0; i_state<n_states; i_state++) {
          for (cell_cnt=0; cell_cnt<global_param.ncells; cell_cnt++)
              free_dist_prcp(&prcp_states[i_state][cell_cnt],veg_con_all[cell_cnt][0].vegetat_type_num);
          free(init_STILL_STORM_states[i_state]);
          free(init_DRY_TIME_states[i_state]);
      }
      free(dmy_states);
  }
  for (cell_cnt=0; cell_cnt<global_param.ncells; cell_cnt++) {    
    free_dist_prcp(&prcp_all[cell_cnt],veg_con_all[cell_cnt][0].vegetat_type_num);
    free_vegcon(&veg_con_all[cell_cnt]);
    free_soilcon(&soil_con_all[cell_cnt]);
  }
  free(veg_con_all);
  free(soil_con_all);
  free(init_STILL_STORM_all);
  free(init_DRY_TIME_all);

  return EXIT_SUCCESS;
}        /* End Main Program */
