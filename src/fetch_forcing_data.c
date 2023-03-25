#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>
 
static char vcid[] = "$Id: fetch_forcing_data.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

double **fetch_forcing_data(float               **nldas_forcing,
                            global_param_struct   global_param,
                            int                   cell_cnt)
/**********************************************************************
  fetch_forcing_data()    Ming Pan, mpan@princeton.edu, May 2005

  This subroutine extracts the forcing for the current grid cell
  from a big forcing array read at the beginning of model run.

**********************************************************************/
{
  extern option_struct    options;
  extern param_set_struct param_set;
  extern int              NR, NF;

  char                 errorstr[MAXSTRING];
  int                  i, i_rec;
  double             **forcing_data;

  /** Allocate data arrays for input forcing data **/
  forcing_data = (double **)calloc(N_FORCING_TYPES,sizeof(double*));
  
  for(i=0;i<N_FORCING_TYPES;i++) 
    if (param_set.TYPE[i].SUPPLIED) 
      forcing_data[i] = (double *)calloc((global_param.nrecs * NF),
			   sizeof(double));

  /** Fetch the current cell from the big forcing array **/
  
  for (i=0; i<N_FORCING_TYPES; i++) {
      
      if (param_set.TYPE[i].SUPPLIED) {
          
          for (i_rec=0; i_rec<global_param.nrecs; i_rec++)
              forcing_data[i][i_rec] = (double) nldas_forcing[i][global_param.nrecs*cell_cnt+i_rec];
      }
  }
  
  return(forcing_data);

}
