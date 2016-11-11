#ifndef READINZGRID_H_INCLUDED
#define READINZGRID_H_INCLUDED



int getNx(FILE *fp);
int getNy(FILE *fp);
void readInZgrid(FILE *fzgrid, double *xgridvals, double *ygridvals, 
				double **xyh);


#endif 
