#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id: make_dist_prcp.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

dist_prcp_struct make_dist_prcp(int  nveg,
				int *Nnodes)
/**********************************************************************
	read_dist_prcp	Keith Cherkauer		May 21, 1996

  This routine creates an array of structures which will store 
  necessary information about the distribution of precipitation, moisture,
  evaporation, and dew.  Mu represents the fractional area of the grid 
  that receives precipitation (wet), while 1-mu is the corresponding 
  area that receives no precipitation.  The value of mu changes with
  the intensity of incoming precipitation, and is set in the routine
  dist_prec.

**********************************************************************/
{
  extern option_struct options;

  dist_prcp_struct temp;
  int              i;

  temp.mu     = (double *)calloc(nveg+1,sizeof(double));
  for ( i = 0; i < nveg + 1; i++ ) temp.mu[i] = 1;
  temp.snow   = make_snow_data(nveg+1);
  temp.energy = make_energy_bal(nveg+1,Nnodes);
  for(i=0;i<2;i++) {
    temp.veg_var[i]  = make_veg_var(nveg);
    temp.cell[i]     = make_cell_data(nveg+1,options.Nlayer);
  }

  return (temp);

}

void copy_dist_prcp(dist_prcp_struct *prcp1,
                    dist_prcp_struct *prcp2,
                    int               nveg) {
    
    extern option_struct options;
    
    int i_dist, i_veg, i_band, i_layer, i_cell;
    
    for (i_veg=0; i_veg<nveg+1; i_veg++) {
        prcp2->mu[i_veg] = prcp1->mu[i_veg];
    }
    
    for (i_veg=0; i_veg<nveg+1; i_veg++) {
        for (i_band=0; i_band<options.SNOW_BAND; i_band++) {
            prcp2->snow[i_veg][i_band]   = prcp1->snow[i_veg][i_band];
            prcp2->energy[i_veg][i_band] = prcp1->energy[i_veg][i_band];
        }
    }
    
    for (i_dist=0; i_dist<2; i_dist++) {
        for (i_veg=0; i_veg<nveg+1; i_veg++) {
            for (i_band=0; i_band<options.SNOW_BAND; i_band++) {
                if (i_veg<nveg) prcp2->veg_var[i_dist][i_veg][i_band] = prcp1->veg_var[i_dist][i_veg][i_band];
                prcp2->cell[i_dist][i_veg][i_band] = prcp1->cell[i_dist][i_veg][i_band];
            }
        }
     }
}

    