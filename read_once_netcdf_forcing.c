#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>
#include <netcdf.h>
 
static char vcid[] = "$Id: read_once_ldas_forcing.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

float **read_once_netcdf_forcing(filenames_struct     *filenames)
/**********************************************************************
  read_once_netcdf_forcing()   Ming Pan, mpan@princeton.edu  May 2017

  This subroutine reads all grid forcing data in one shot.

**********************************************************************/
{
  extern option_struct    options;
  extern param_set_struct param_set;
  extern int              NR, NF;
  extern int              flag;
  
  extern global_param_struct   global_param;
  
  FILE         *fp_soil;
  char         fnin[MAXSTRING], dmychar[MAXSTRING], tmpstr[MAXSTRING];
  int          dummy, ncells, nvars;
  int          grid_i, grid_j, i_cell, i_rec;
  float        *tmp_float;
  double       *tmp_double;

  char                 errorstr[MAXSTRING];
  int                  i;
  long                 nread;
  float             **forcing_data;
  
  int          ncid, varid, retval, i_var;
  char         ncvarnames[N_FORCING_TYPES][MAXSTRING], ncfilenames[N_FORCING_TYPES][MAXSTRING];
  int          provided[N_FORCING_TYPES];
  size_t       start[3], count[3];
  
  /* Allocate memory for tmp_grid */
  if (param_set.FORCE_FORMAT[0] == NCDOUBLE) {
      tmp_double = (double *) calloc(options.INPUT_GRID_NX*options.INPUT_GRID_NY, sizeof(double));
      if (tmp_double == NULL)
          vicerror("Memory allocation error in read_once_netcdf_forcing(): tmp_double");
  } else {
      tmp_float = (float *) calloc(options.INPUT_GRID_NX*options.INPUT_GRID_NY, sizeof(float));
      if (tmp_float == NULL)
          vicerror("Memory allocation error in read_once_netcdf_forcing(): tmp_float");
  }
  
  /* Read soil param file to find out # of cells and lat-lons */
  
  fp_soil = open_file(filenames->soil, "r");
      /*vicerror("Can't open soil parameter file read_once_netcdf_forcing().");*/
  
  fgets(tmpstr, MAXSTRING, fp_soil);
  
  ncells = 0;
  
  while (!feof(fp_soil)) {
      
      sscanf(tmpstr, "%d %d %lf %lf", &flag, &dummy, &global_param.lat[ncells], &global_param.lon[ncells]);
      
      if (flag) ncells++;
      
      fgets(tmpstr, MAXSTRING, fp_soil);
      
  }
  
  fclose(fp_soil);
  global_param.ncells = ncells;
  
  nvars = 0;
  
  /** Allocate data arrays for input forcing data **/
  forcing_data = (float **)calloc(N_FORCING_TYPES,sizeof(float*));
  for(i=0;i<N_FORCING_TYPES;i++) {
      if (param_set.TYPE[i].SUPPLIED) {
          provided[nvars] = i;
          strcpy(ncvarnames[nvars], param_set.TYPE[i].ncvarname);
          strcpy(ncfilenames[nvars], param_set.TYPE[i].ncfilename);
          forcing_data[i] = (float *)calloc((global_param.nrecs * ncells), sizeof(float));
          if (forcing_data[i] == NULL)
              vicerror("Memory allocation error in read_once_netcdf_forcing(): forcing_data");
          nvars++;
      }
  }
  
  fprintf(stderr, "Forcing: Ncells = %d, Nvars = %d, Nrecs = %d\nTotal # of floats/bytes allocated: %d/%ld\n",
              ncells, nvars, global_param.nrecs, ncells*nvars*global_param.nrecs, ncells*nvars*global_param.nrecs*sizeof(float));

  /* read files */
  for (i_var=0; i_var<nvars; i_var++) {
      
      fprintf(stderr, "Reading forcing variable #%d ...\n", i_var+1);
      sprintf(dmychar, "_%4d%02d", global_param.forceyear[0], global_param.forcemonth[0]);
      strcpy(fnin, filenames->forcing[0]);
      strcat(fnin, "/");
      strcat(fnin, ncfilenames[i_var]);
      //strcat(fnin, dmychar);
      strcat(fnin, ".nc");
      
      if ((retval = nc_open(fnin, NC_NOWRITE, &ncid))) {
          sprintf(errorstr, "Cannot open forcing file %s.", fnin);
          vicerror(errorstr);
      }
      if ((retval = nc_inq_varid(ncid, ncvarnames[i_var], &varid))) {
          sprintf(errorstr, "Cannot find forcing variable %s.", ncvarnames[i_var]);
          vicerror(errorstr);
      }
      
      start[0] = 0; start[1] = 0; start[2] = 0;
      count[0] = 1; count[1] = options.INPUT_GRID_NY; count[2] = options.INPUT_GRID_NX;
      
      fprintf(stderr, "Skipping %d records.\n", global_param.forceskip[0]);
      for (i_rec=0; i_rec<global_param.nrecs; i_rec++) {
          
          start[0] = i_rec + global_param.forceskip[0];
                    
          if (param_set.FORCE_FORMAT[0] == NCDOUBLE) {
              if ((retval = nc_get_vara_double(ncid, varid, start, count, tmp_double))) {
                  sprintf(errorstr, "Forcing file reading error for variable %s.", ncvarnames[i_var]);
                  vicerror(errorstr);
              }
          } else {
              if ((retval = nc_get_vara_float(ncid, varid, start, count, tmp_float))) {
                  sprintf(errorstr, "Forcing file reading error for variable %s.", ncvarnames[i_var]);
                  vicerror(errorstr);
              }
          }
          
          for (i_cell=0; i_cell<ncells; i_cell++) {
              
              grid_i = (global_param.lon[i_cell]-options.INPUT_GRID_XO)/options.INPUT_GRID_DX + 0.5;
              grid_j = (global_param.lat[i_cell]-options.INPUT_GRID_YO)/options.INPUT_GRID_DY + 0.5;
              
              if (options.YREV) grid_j = options.INPUT_GRID_NY-1 - grid_j;
              
              if (param_set.FORCE_FORMAT[0] == NCDOUBLE)
                  forcing_data[provided[i_var]][global_param.nrecs*i_cell+i_rec] = (float) tmp_double[grid_j*options.INPUT_GRID_NX+grid_i] * param_set.TYPE[provided[i_var]].multiplier;
              else
                  forcing_data[provided[i_var]][global_param.nrecs*i_cell+i_rec] = tmp_float[grid_j*options.INPUT_GRID_NX+grid_i] * param_set.TYPE[provided[i_var]].multiplier;
              
          }
          
      } /* end of rec loop */
      
      if ((retval = nc_close(ncid))) {
          sprintf(errorstr, "Cannot close forcing file %s.", fnin);
          vicerror(errorstr);
      }
          
  } /* end of var loop */
      
  if (param_set.FORCE_FORMAT[0] == NCDOUBLE) free(tmp_double);
  else free(tmp_float);

  return(forcing_data);

}
