#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
#include <netcdf.h>

static char vcid[] = "$Id: write_data_nldas_grid.c,v 4.1.2.4 2004/05/06 20:26:41 tbohn Exp $";

void write_data_netcdf(out_data_struct     *out_data,
                       filenames_struct    *filenames,
                       global_param_struct *global,
                       int                  cell_cnt,
                       int                  nrecs,
                       int                  rec)
/**********************************************************************
  write_data_netcdf(), Ming Pan, May 2017, modified from write_data()

  This subroutine writes all energy and moisture balance parameters to
  a single gridded output file.

  OUTPUT:
	evaporation and vapor fluxes in mm/time step
	layer moisture in mm/time step
	runoff in mm/time step
	baseflow in mm/time step
	freezing and thawing depths in cm
	snow depth in cm
	snow water equivlence in mm
	all energy fluxes are in W/m^2

  Modifications:
  5/20/96	Program was modified to account for a variable
		number of soil layers.  It was also modified to
		write out frozen soils data per time step.	KAC
  1/15/97	Program modified to output daily sums, or values
		independant of selected time step.  This aids in
		comparisons between model versions.		KAC
  3/98          Routine modified to output fluxes in PILPS2c 
                ASCII column format                             Dag
  4/30/98       Routine modified to add binary output options for
                improved file speed, and less disk usage for large
		model basins                                    KAC
  7/19/99       modified to output a single binary file containing
                the data selected for the LDAS project         KAC
  8/3/99        modified again to reduce the storage space needed
                for the LDAS output files.  
  1/4/2000      modified to allow both standard and LDAS formatted
                output using a compiler flag                    KAC
  3-12-03   added energy fluxes to snow band output files   KAC
  04-23-2003    modified LDAS SWQ output, so that it is multiplied by
                10 instead of 100 before being converted to a short
                integer.  This reduces stored value precision to 0.1,
                but increases the maximum storable SWQ, which was
                exceeded in previous LDAS simulations.          KAC

**********************************************************************/
{
  extern option_struct options;

  static int          nvars, ncells;
  static float       *tmp_grid;
  static double       *lats, *lons;
  
  int                 i, j, k;
  int                 grid_i, grid_j, i_var, i_cell, i_rec;
  char                tmpstr[MAXSTRING];
  int                 dummy;
  float               *tmp_out, smtot, runtot;
  
  char                fnout[MAXSTRING], dmychar[MAXSTRING], monthchar[MAXSTRING], timechar[MAXSTRING], errorstr[MAXSTRING];
  
  int                 tmpyear, tmpmonth, tmpday, tmphour, tmpjday=0;
  long                output_bytes;

  int     retval, ncid, lon_dimid, lat_dimid, rec_dimid;
  int     lat_varid, lon_varid, rec_varid, varid;
  int     dimids[3];
  size_t  start[3], count[3];
  int     *timesteps;
  static float undef[] = {UNDEF};
  
  char   varnames[50][MAXSTRING], varlong[50][MAXSTRING], varunits[50][MAXSTRING];
  char  *global_model = "VIC 4.0.5 daily image parallel version, re-coded by Ming Pan, mpan@princeton.edu, fallspinach@gmail.com";
  /*
  char  *global_project = "EDgE";
  char  *global_producer = "Ming Pan <mpan@princeton.edu>";
  */
  char  *global_title = "VIC model simulations";
  
  //fprintf(stderr, "Start writing cell %d, record %d\n", cell_cnt, rec);
  
  /* initialization */
  
  if (rec == -1) {
      
      ncells = global->ncells;
      
      /* Determine # of variables to output */
      
      nvars  = 2;                                                  /* evap, runoff+baseflow */
      nvars += 1;                                                  /* column total moist */
      nvars += 1;                                                  /* swq */
      
      strcpy(varnames[0], "evap");     strcpy(varunits[0], "m/day"); strcpy(varlong[0], "Total Evapotranpiration");
      strcpy(varnames[1], "runoff");   strcpy(varunits[1], "m/day"); strcpy(varlong[1], "Total Runoff");
      strcpy(varnames[2], "sm");       strcpy(varunits[2], "m");     strcpy(varlong[2], "Column Total Soil Moisture");
      strcpy(varnames[3], "swe");      strcpy(varunits[3], "m");     strcpy(varlong[3], "Snow Water Equivalent");
      
      /* Allocate memory */
      
      fprintf(stderr, "Output: Ncells = %d, Nvars = %d, Nrecs = %d\nTotal # of floats/bytes to be allocated: %d/%ld\n",
              ncells, nvars, nrecs, ncells*nvars*nrecs, ncells*nvars*nrecs*sizeof(float));
      
      tmp_grid = (float *) calloc(ncells*nvars*nrecs, sizeof(float));
      if (tmp_grid == NULL)
          vicerror("Memory allocation error in write_data_netcdf().");
      
      return;
  }
  
  /* write to file */
  
  if (rec == -2) {
      
      /** Check output file size **/
      output_bytes = options.OUTPUT_GRID_NX*options.OUTPUT_GRID_NY*nrecs*sizeof(float);
      fprintf(stderr, "Output data grid size %ld bytes per variable.\n", output_bytes);
      
      /* Allocate memory for tmp_out */
      
      tmp_out = (float *) calloc(options.OUTPUT_GRID_NX*options.OUTPUT_GRID_NY, sizeof(float));
      if (tmp_out == NULL)
          vicerror("Memory allocation error in write_data_netcdf(): tmp_out.");
      
      timesteps = (int *) calloc(nrecs, sizeof(int));
      if (timesteps == NULL)
          vicerror("Memory allocation error in write_data_netcdf(): timesteps.");
      for (i_rec=0; i_rec<nrecs; i_rec++) timesteps[i_rec] = i_rec;
      
      lats = (double *) calloc(options.OUTPUT_GRID_NY, sizeof(double));
      lons = (double *) calloc(options.OUTPUT_GRID_NX, sizeof(double));
      for (i=0; i<options.OUTPUT_GRID_NY; i++) lats[i] = options.OUTPUT_GRID_YO + i*options.OUTPUT_GRID_DY;
      for (j=0; j<options.OUTPUT_GRID_NX; j++) lons[j] = options.OUTPUT_GRID_XO + j*options.OUTPUT_GRID_DX;

      /* write output rec by rec, var by var */
      
      fprintf(stderr, "\nStart writing to file.\n");
      
      tmpyear  = global->startyear;
      tmpmonth = global->startmonth;
      tmpday   = global->startday;
      tmphour  = global->starthour;

      /* open output file */
      sprintf(dmychar, "_%04d%02d", tmpyear, tmpmonth);
      sprintf(timechar, "days since %04d-%02d-%02d 00:00", tmpyear, tmpmonth, tmpday);
      
      for (i_var=0; i_var<nvars; i_var++) {
          
          strcpy(fnout, filenames->result_dir);
          strcat(fnout, varnames[i_var]);
          strcat(fnout, "_VIC");
          strcat(fnout, dmychar);
          strcat(fnout, ".nc");
          
          /* define dimensions */
          if ((retval = nc_create(fnout, NC_CLOBBER|NC_NETCDF4, &ncid))) {
              sprintf(errorstr, "Cannot create output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_def_dim(ncid, "lat", options.OUTPUT_GRID_NY, &lat_dimid))) {
              sprintf(errorstr, "Cannot create lat dimension in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_def_dim(ncid, "lon", options.OUTPUT_GRID_NX, &lon_dimid))) {
              sprintf(errorstr, "Cannot create lon dimension in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_def_dim(ncid, "time", NC_UNLIMITED, &rec_dimid))) {
              sprintf(errorstr, "Cannot create time (record) dimension in output file %s.", fnout);
              vicerror(errorstr);
          }
          
          /* define dimension variables */
          if ((retval = nc_def_var(ncid, "lat", NC_DOUBLE, 1, &lat_dimid, &lat_varid))) {
              sprintf(errorstr, "Cannot create lat variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_def_var(ncid, "lon", NC_DOUBLE, 1, &lon_dimid, &lon_varid))) {
              sprintf(errorstr, "Cannot create lon variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_def_var(ncid, "time", NC_INT, 1, &rec_dimid, &rec_varid))) {
              sprintf(errorstr, "Cannot create time variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_put_att_text(ncid, lat_varid, "units", strlen("degrees_north"), "degrees_north"))) {
              sprintf(errorstr, "Cannot define units for lat variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_put_att_text(ncid, lon_varid, "units", strlen("degrees_east"), "degrees_east"))) {
              sprintf(errorstr, "Cannot define units for lon variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_put_att_text(ncid, rec_varid, "units", strlen(timechar), timechar))) {
              sprintf(errorstr, "Cannot define units for time variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_put_att_text(ncid, rec_varid, "calendar", strlen("standard"), "standard"))) {
              sprintf(errorstr, "Cannot define calendar for time variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_put_att_text(ncid, lat_varid, "long_name", strlen("latitude"), "latitude"))) {
              sprintf(errorstr, "Cannot define long_name for lat variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_put_att_text(ncid, lon_varid, "long_name", strlen("longitude"), "longitude"))) {
              sprintf(errorstr, "Cannot define long_name for lon variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_put_att_text(ncid, rec_varid, "long_name", strlen("time"), "time"))) {
              sprintf(errorstr, "Cannot define long_name for time variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          
          /* define data varialbes */
          dimids[0] = rec_dimid;
          dimids[1] = lat_dimid;
          dimids[2] = lon_dimid;
          
          if ((retval = nc_def_var(ncid, varnames[i_var], NC_FLOAT, 3, dimids, &varid))) {
              sprintf(errorstr, "Cannot create data variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_put_att_text(ncid, varid, "units", strlen(varunits[i_var]), varunits[i_var]))) {
              sprintf(errorstr, "Cannot define units for data variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_put_att_text(ncid, varid, "long_name", strlen(varlong[i_var]), varlong[i_var]))) {
              sprintf(errorstr, "Cannot define long_name for data variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_put_att_float(ncid, varid, "_FillValue", NC_FLOAT, 1, undef))) {
              sprintf(errorstr, "Cannot define _FillValue for data variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          /* global attributes */
          if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "model", strlen(global_model), global_model))) {
              sprintf(errorstr, "Cannot define global attribute 'model' in output file %s.", fnout);
              vicerror(errorstr);
          }
          /*
          if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "project", strlen(global_project), global_project))) {
              sprintf(errorstr, "Cannot define global attribute 'project' in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "producer", strlen(global_producer), global_producer))) {
              sprintf(errorstr, "Cannot define global attribute 'producer' in output file %s.", fnout);
              vicerror(errorstr);
          }
          */
          if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "title", strlen(global_title), global_title))) {
              sprintf(errorstr, "Cannot define global attribute 'title' in output file %s.", fnout);
              vicerror(errorstr);
          }
          
          if ((retval = nc_def_var_deflate(ncid, varid, 0, 1, 5))) {
              sprintf(errorstr, "Cannot deflate data variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          
          if ((retval = nc_enddef(ncid))) {
              sprintf(errorstr, "Cannot end def mode in output file %s.", fnout);
              vicerror(errorstr);
          }
          
          /* write variables */
          if ((retval = nc_put_var_double(ncid, lat_varid, lats))) {
              sprintf(errorstr, "Cannot write lat variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          if ((retval = nc_put_var_double(ncid, lon_varid, lons))) {
              sprintf(errorstr, "Cannot write lon variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          
          start[0] = 0; start[1] = 0; start[2] = 0;
          count[0] = 1; count[1] = options.INPUT_GRID_NY; count[2] = options.INPUT_GRID_NX;
          
          for (i_rec=0; i_rec<nrecs; i_rec++) {
              
              start[0] = i_rec;
                            
              /* initialize grid */
              for (i=0; i<options.OUTPUT_GRID_NY; i++)
                  for (j=0; j<options.OUTPUT_GRID_NX; j++) tmp_out[options.OUTPUT_GRID_NX*i+j] = UNDEF;
              
              for (i_cell=0; i_cell<ncells; i_cell++) {
                  
                  /* find out where the cell is */
                  grid_i = (global->lon[i_cell]-options.OUTPUT_GRID_XO)/options.OUTPUT_GRID_DX + 0.5;
                  grid_j = (global->lat[i_cell]-options.OUTPUT_GRID_YO)/options.OUTPUT_GRID_DY + 0.5;
                  
                  tmp_out[options.OUTPUT_GRID_NX*grid_j+grid_i] = tmp_grid[nvars*ncells*i_rec+ncells*i_var+i_cell];
                  
              }
              
              if ((retval = nc_put_vara_float(ncid, varid, start, count, tmp_out))) {
                  sprintf(errorstr, "Cannot write data variable in output file %s.", fnout);
                  vicerror(errorstr);
              }
              
          } // end of rec loop
          
          if ((retval = nc_put_var_int(ncid, rec_varid, timesteps))) {
              sprintf(errorstr, "Cannot write time variable in output file %s.", fnout);
              vicerror(errorstr);
          }
          
          if ((retval = nc_close(ncid))) {
              sprintf(errorstr, "Cannot close output file %s.", fnout);
              vicerror(errorstr);
          }
                    
      } // end of var loop
      
      free(tmp_grid);
      free(timesteps);
      free(lats); free(lons);
      
      return;
  }
  
  /************************************
    Attention: NO Frozen Soil or Snow Band Variables will be in the output
  ************************************/
  
  /************************************
    Output Standard Energy and Moisture Flux Variables
  ************************************/
  
  i_var = 0;
  
  //fprintf(stderr, "started writing at byte: %d ", nvars*ncells*rec+ncells*i_var+cell_cnt);
  
  /***** Write Binary Fluxes Variables *****/
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt] = (float) (out_data->evap/1000.0);      i_var++;
  runtot = out_data->runoff + out_data->baseflow;
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt] = (float) (runtot/1000.0);              i_var++;
  smtot = 0;
  for(j=0;j<options.Nlayer;j++) {    
    smtot += (float)out_data->moist[j];
  }
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt] = (float) (smtot/1000.0);               i_var++;
  /***** Write Binary Snow Variables *****/
  tmp_grid[nvars*ncells*rec+ncells*i_var+cell_cnt] = (float) (out_data->swq[0]/1000.0);    i_var++;
  //fprintf(stderr, "finish writing cell %d, record %d\n", cell_cnt, rec);

}
