#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "matrixOperations.h"
#include "f2c.h"


#define MAXNAME 100
#define MAXLINEBUFFER 5000

int getNx(FILE *fp){
	char linebuffer[MAXLINEBUFFER];
	int nx;
	fgets(linebuffer, MAXLINEBUFFER, fp);
	int c = sscanf(linebuffer, "Nx=%d", &nx);
	//int c = fscanf(fp, "Nx=%d", &nx);
	//printf("%d item(s) matched\n", c);
	rewind(fp);
	return nx;
} 

int getNy(FILE *fp){
	char linebuffer[MAXLINEBUFFER];
	int ny;
	fgets(linebuffer, MAXLINEBUFFER, fp); // read 1st line
	fgets(linebuffer, MAXLINEBUFFER, fp); // read 2nd line
	int c = sscanf(linebuffer, "Ny=%d", &ny);
	//int c = fscanf(fp, "Nx=%d", &nx);
	//printf("%d item(s) matched\n", c);
	rewind(fp);
	return ny;
} 




void readInZgrid(FILE *fzgrid, double *xgridvals, double *ygridvals, 
				double **xyh){
	
		
	
	char linebuffer[MAXLINEBUFFER];

	int Nx, Ny;
	double x1, x2, y1, y2;
	
	int i;
	char *pos;
	int counter;
	int *nxny = (int*) malloc(sizeof(int) * 2);
	double *xyparms = (double*) malloc(sizeof(double) * 4);
	int tmpint;
	double tmpdouble;
	int xyparmsInd = 0;
	
	for (i = 0; i < 2; i++){
		counter = 0;
		fgets(linebuffer, MAXLINEBUFFER, fzgrid);
		//puts(linebuffer);
		pos = strtok(linebuffer," ");
		
	
		while (pos != NULL){
			//puts(pos);
			if (counter == 0){
				if (i == 0){
					sscanf(pos, "Nx=%d:", &tmpint);
				}
				else{
					sscanf(pos, "Ny=%d:", &tmpint);
				}
				nxny[i] = tmpint;
			}
		
			if (counter == 2){
				sscanf(pos, "%lf,", &tmpdouble);
				//printf("x1 = %7.3f\n", tmpdouble);
				xyparms[xyparmsInd] = tmpdouble;
				xyparmsInd++;
			}
			if (counter == 3){
				sscanf(pos, "%lf,", &tmpdouble);
				//printf("x1 = %7.3f\n", tmpdouble);
				xyparms[xyparmsInd] = tmpdouble;
				xyparmsInd++;
			}
			pos = strtok(NULL, " ");
			counter++;
			//printf("-----------------\n");
		}
		
	} 
	//vect_Print(xyparms, 4, 4);
	//vect_PrintINT(nxny, 2);
	int nx = nxny[0]; int ny = nxny[1]; free(nxny);
	printf("nx = %d\t ny = %d\n", nx, ny);
	x1 = xyparms[0]; x2 = xyparms[1]; y1 = xyparms[2]; y2 = xyparms[3];
	free(xyparms);
	
	/* Remove either 'Elevation=' or 'Pileheight=' */
	fgets(linebuffer, MAXLINEBUFFER, fzgrid);
	//puts(linebuffer);
	
	
	
	int j;
	for (i = 0; i<ny; i++){
		fgets(linebuffer, MAXLINEBUFFER, fzgrid);
		pos = strtok(linebuffer, " ");
		for(j =0; j<nx; j++){
			xyh[i][j] = atof(pos);
			//printf("(%d, %d) \t %3.9f\n", i,j,xyh[i][j]);
			pos = strtok(NULL, " ");
		}
	}
	
	for(i= 0; i<nx; i++){
		xgridvals[i] = ((((2*i) + 0.5)/(2*nx)) * (x2-x1)) + x1;
	}
	
	for(i = 0; i<ny; i++){
		ygridvals[i] = ((((2*i) + 0.5)/(2*ny)) * (y2-y1)) + y1;
	
	}	

}


