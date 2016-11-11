#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>



#include "matrixOperations.h"
#include "readInZgrid.h"

#define MAXNAME 100
#define MAXLINEBUFFER 5000


int linesInFile(FILE *fp){	

    char buffer[MAXLINEBUFFER];
    int nrows = 0;
    while (fgets(buffer, MAXLINEBUFFER, fp)){
        nrows++;
    }
    rewind(fp);
	return nrows;
}

int colsInFile(FILE *fp){
	
	char buffer[MAXLINEBUFFER];
	char *pos;
	int counter = 0;
	fgets(buffer, MAXLINEBUFFER, fp);
	pos = strtok(buffer, " ");
	
	while (pos != NULL){
		pos = strtok(NULL, " ");
		counter++;
	}
	rewind(fp);
	return counter;
}

void readVectFromFileINT (FILE *fp, int *vect, int size){
	
	int i;
	for (i = 0; i< size; i++){
		fscanf(fp, "%d", &vect[i]);
	}
}

int main(int argc, char *argv[]){
	
	clock_t tic = clock();
	
	/*  Read in indlistnorth and indlistwest */
	char fnameIndLstN[MAXNAME];
	char fnameIndLstW[MAXNAME];
	strcpy(fnameIndLstN, "./newrunsdir_BV_100m/indlistnorth.txt");
	strcpy(fnameIndLstW, "./newrunsdir_BV_100m/indlistwest.txt");
    puts(fnameIndLstW);
    FILE *findlstN;
	FILE *findlstW;
	findlstN = fopen(fnameIndLstN, "r");
	findlstW = fopen(fnameIndLstW, "r");
	int nbrIndLst = linesInFile(findlstN);
	//printf("%d inds\n", nbrIndLst);
	
	int *indN = (int*) malloc(sizeof(int) * nbrIndLst);
	int *indW = (int*) malloc(sizeof(int) * nbrIndLst);
	readVectFromFileINT(findlstN, indN, nbrIndLst);
	readVectFromFileINT(findlstW, indW, nbrIndLst);
	fclose(findlstN);
	fclose(findlstW);
	//vect_PrintINT(indW, nbrIndLst);
	
	char dirname[MAXNAME];
	strcpy(dirname,"./oct_output_files_20cells"); 
	char filename[MAXNAME];
	int irun;
	
	// file-dependent variables
	FILE *fzgrid;
	double **xyh;
	double *xgridvals, *ygridvals;
	double dxdy, vol, minmaxh;
	
	
	/* read elevation_20 grid*/
	FILE *felev20;
	felev20 = fopen("./newrunsdir/elevation_20.grid", "r");
	int elevNx = getNx(felev20);
	int elevNy = getNy(felev20);
	//printf("nx = %d , ny = %d \n", elevNx, elevNy);
	double **elev20grid;
	elev20grid = create2Darray(elevNy, elevNx, 0.0);
	double *elevxgrid = (double*) malloc(sizeof(double) * elevNx);
	double *elevygrid = (double*) malloc(sizeof(double) * elevNy);
	readInZgrid(felev20, elevxgrid, elevygrid, elev20grid); 
	

	int nbrRuns = 2048;
	
	int i, j;
	double msum;
	double **heights= create2Darray(elevNx, elevNy, 0.0);
    double *hsaveCol = (double*) calloc(nbrIndLst,sizeof(double));
	

    /*BEGINNING OF ITERATIONS */
    irun = atoi(argv[1]) + 2 ; // pileheight* start at 2	
	
    snprintf(filename, MAXNAME, "%s/pileheightrecord.%06d", dirname,irun);
    puts(filename);

    fzgrid = fopen(filename, "r");
    int nx = getNx(fzgrid);
    int ny = getNy(fzgrid);
    //printf("nx = %d , ny = %d \n", nx, ny);
    xyh = create2Darray(ny, nx, 0.0);
    xgridvals = (double*) malloc(sizeof(double) * nx);
    ygridvals = (double*) malloc(sizeof(double) * ny);
    readInZgrid(fzgrid, xgridvals, ygridvals, xyh);
    fclose(fzgrid);

    dxdy = (xgridvals[1] - xgridvals[0])*(ygridvals[1]- ygridvals[0]);
    //printf("dxdy = %5.6f\n", dxdy);
    
    vol = matrix_Sum(xyh, ny, nx) * dxdy;
    //printf("vol = %e\n", vol);
    
    vect2DInterp2(xgridvals, ygridvals, xyh, elevxgrid, elevygrid,
                    heights, nx, ny, elevNx, elevNy);
    
    j = 0;
    for (i = 1; i< elevNy; i+= 2){
        heights[elevNx-1][i] = xyh[j][nx-1];
        heights[0][i] = xyh[j][0];
        j++;
            
    }

    for (i = 2; i<elevNy; i+= 2){
        heights[elevNx-1][i] = (heights[elevNx-1][i-1] + 
                                heights[elevNx-1][i+1])/2;
        heights[0][i-2] = (heights[elevNx-1][i-1] + 
                            heights[elevNx-1][i+1])/2;
    }
    
    j = 0;
    for(i = 1; i<elevNx; i+=2){
        heights[i][elevNy-1] = xyh[ny-1][j];
        heights[i][0] = xyh[0][j];
        j++;
    }
    
    for(i= 2; i<elevNx; i+= 2){
        heights[i][elevNy-1] = (heights[i-1][elevNy-1] + 
                                heights[i+1][elevNy-1])/2;
        heights[i][0] = (heights[i-1][elevNy-1] + 
                        heights[i+1][elevNy-1])/2;
    }
    
    for(j = 0; j<nbrIndLst; j++){ 
        hsaveCol[j] = heights[indW[j]-1][indN[j]-1];
    }
    

    /* END OF FOR IRUN LOOP */ 


	printf("***** \t \t  Finished run: %f secs \t \t ***** \n\n",
		   (double)(clock() - tic) / CLOCKS_PER_SEC);

    strcpy(dirname, "./newrunsdir");
    snprintf(filename, MAXNAME, "%s/heightsCols_%06d.txt", dirname, irun);
    puts(filename);
    FILE *fhcol;
    fhcol = fopen(filename, "w");
    vect_PrintToFileDOUBLE(hsaveCol, fhcol, nbrIndLst,6);
    fclose(fhcol);


	return 0;
}
