#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>



/* LAPACK & MYINCLUDES  */
#include "matrixOperations.h"
#include "prTimesLhood.h"
#include "f2c.h" // NB: main.c only needs this for abs, min, and max.


#define MAXLINEBUFFER 200
#define MAXNAME 200
//void MAIN__(){}
//void MAIN_(){}
//void _MAIN_(){}



int linesInFile(FILE *fp){	


    char buffer[MAXLINEBUFFER];
    int nrows = 0;
    while (fgets(buffer, MAXLINEBUFFER, fp)){
        nrows++;
    	//printf("(%d)\t %s", nrows, buffer);

    }
    //printf("nrows= %d\n", nrows);
    rewind(fp);
	return nrows;
}






void readDesignPoints(double* vLst, double* ang0Lst,
                      double* basFrLst, double* intAngLst,
                      double* volTimesBasLst, int npts, FILE* finput){
    int maxcols = 200;
    char line[maxcols];
    char* pos;
    int counter;
    int rownbr = 0;
    while (fgets(line, maxcols, finput) != NULL){
        //puts(line);
        pos = strtok(line, " ");
        counter = 0;
        while (pos != NULL){
            //printf("%d ", counter);
            if (counter == 0){
                vLst[rownbr] = atof(pos);
                //printf("vLst");
            }
            else if (counter == 1){
                ang0Lst[rownbr]= atof(pos);
                //printf(" ang0Lst");
            }
            else if (counter == 2){
                basFrLst[rownbr] = atof(pos);
                //printf(" basFrList");
            }
            else if (counter == 3){
                intAngLst[rownbr] = atof(pos);
                //printf(" intAngLst");
            }

            else if (counter == 4){
                volTimesBasLst[rownbr] = atof(pos);
                //printf(" volTimesBasLst");
            }

            //printf("%5.5f ", atof(pos));
            pos = strtok(NULL, " ");
            counter++;
        }
        rownbr++;
        //printf("\n");
        //printf("%d\n", rows);
    }
    //printf("There are %d rows\n", rownbr); // confirmation message.

}

 
void myPrintToFile(double* mat,FILE* foutput,int nrows, int ncols){
    int i;

    int totsize = nrows * ncols;
    printf("Printing to file ...\n");
    for(i = 0; i< totsize; i++){
            if (((i+1)%ncols) == 0){
                fprintf(foutput, "%3.4f\n", *(mat+i));
            }
            else{
                fprintf(foutput, "%3.4f ", *(mat+i));
            }
    }
    printf("Done Printing to file\n");
}

void readLineSamples(FILE *fLines, int nlines, double *linesamples1, 
					double *linesamples2){
	
	
	char buffer[MAXLINEBUFFER]; char *pos; int counter;
    double r1, r2;
	int i;
	int rownbr = 0;
	for(i = 0; i<nlines; i++){
		fgets(buffer, MAXLINEBUFFER, fLines);
		pos = strtok(buffer, " ");
		counter = 0;
		while (pos != NULL){
			if(counter == 0){
				linesamples1[rownbr] = atof(pos);
			}
			else if (counter== 1){
				linesamples2[rownbr] = atof(pos);
			}
			pos = strtok(NULL, " ");
			counter++;
		}
		rownbr++;
	}
    
    //for(i = 0; i<nlines; i++){
    	//printf("%d\t= \t %3.8f \t %3.8f\n",i, linesamples1[i], 
    	//		linesamples2[i]);
   // }
    
    //printf("nlines = %d\n", nlines); 					
}

void readDaniloQs(FILE *fDaniloqs, int nlines, double **daniloQsArray){
	
	char buffer[MAXLINEBUFFER];
	char *pos;
	int colCounter = 0;
	int rowCounter = 0;
	while (fgets(buffer, MAXLINEBUFFER, fDaniloqs) != NULL){
		pos = strtok(buffer, " ");
		colCounter = 0;
		while (pos != NULL){
			daniloQsArray[rowCounter][colCounter] = atof(pos);			
			pos = strtok(NULL, " ");
			colCounter++;
		}
		rowCounter++;
	}

}

//void readQuants(FILE *fq, int nlines, double **quantsArray){
void readQuants(FILE *fq,int nlines, double *vols, double *q1, double *q5,
				double *q50, double *q95, double *q99){
	
	
	char buffer[MAXLINEBUFFER];
	char *pos;
	int colCounter = 0;
	int rowCounter = 0;
	while (fgets(buffer, MAXLINEBUFFER, fq) != NULL){
		pos = strtok(buffer, " ");
		colCounter = 0;
		while (pos != NULL){
			//quantsArray[rowCounter][colCounter] = atof(pos);			
			
			if (colCounter == 0){
				vols[rowCounter] = atof(pos);
			}
			else if (colCounter == 1){
				q1[rowCounter] = atof(pos);
			}
			else if (colCounter == 2){
				q5[rowCounter] = atof(pos);
			}
			else if (colCounter == 3){
				q50[rowCounter] = atof(pos);
			}
			else if (colCounter == 4){
				q95[rowCounter] = atof(pos);
			}
			else if (colCounter == 5){
				q99[rowCounter] = atof(pos);
			}
			pos = strtok(NULL, " ");
			colCounter++;
		}
		rowCounter++;
	}

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



void readVectFromFile(FILE *fp, double *vect, int size){
	
	int i;
	for (i = 0; i< size; i++){
		fscanf(fp, "%lf", &vect[i]);
	}
}



int main(int argc, char *argv[]){
	srand(time(NULL));
	clock_t tic = clock();
	int i, j;
	
	
    /* Read in design points, each column in a separate array*/
    FILE *fxp;
    int idfile = atoi(argv[1]);
    char indirname[MAXNAME];
    strcpy(indirname, "./subdesignsIO");
    char designfilename[MAXNAME];
    snprintf(designfilename, MAXNAME, "%s/designp_%06d.txt", indirname,
            idfile);
    puts(designfilename);
    
    fxp = fopen(designfilename,"r");
    int npts = linesInFile(fxp); // number of design points
    int xcols = colsInFile(fxp);

    if (npts == 1){
        printf("File has 1 design point.\n");
        exit(0);
    }


    printf("Input file has %d rows and %d column(s).\n", npts, xcols);
    double *volumeVect = (double*) malloc(sizeof(double) * npts); 
    double *ang0Vect = (double*) malloc(sizeof(double) * npts); 
    double *basFricAngVect = (double*) malloc(sizeof(double) * npts); 
    double *intAngVect = (double*) malloc(sizeof(double) * npts); 
	double *volBasVect = (double*) malloc(sizeof(double) * npts);
    readDesignPoints(volumeVect, ang0Vect, basFricAngVect,
                    intAngVect, volBasVect, npts, fxp);
    fclose(fxp);
    

    
    
    /* Read in logheights values*/
    FILE *fheights;
    char logheightfilename[MAXNAME];
    snprintf(logheightfilename, MAXNAME, "%s/logheightp_%06d.txt",
                indirname, idfile);
    fheights=  fopen(logheightfilename, "r");
    double *y = (double*) malloc(sizeof(double) * npts);
    readVectFromFile(fheights, y, npts);
    fclose(fheights);
	//vect_Print(y, npts, 4);

    /* Read in bofv_v*/
    FILE *fbofv;
    fbofv = fopen("bofv_v.txt", "r");
    int bofv_len = linesInFile(fbofv);
    double *bofv_v = (double*) malloc(sizeof(double) * bofv_len); 
    readVectFromFile(fbofv, bofv_v, bofv_len);
    fclose(fbofv);

    
    
    double *pars = (double*) malloc(sizeof(double)* 2);
    vect_LinFit(volumeVect, y, npts, pars);
    //vect_Print(pars, 2, 4);

    double m_v = pars[1]; // slope term
    double m_bf = 0;
    double m_vbf = 0;
    double beta_v = pars[0]; // y-intercept term
    double parest[] = {beta_v, m_v, 0, 0};
    double circterms[] = {1,2,3,4,5};

    int ndp = npts;
    
    int Npred = 100; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    //converting angles from radians to degrees
    vect_MultN(ang0Vect, ang0Vect, PI/180, npts);
    double *vgrid = (double*) malloc(sizeof(double) * Npred);
    double atemp = min(vect_Max(volumeVect, ndp)-0.5, 5.02);
    double btemp = min(vect_Max(volumeVect, ndp)+0.25, 9.497);
    vect_Linspace(vgrid, atemp, btemp, Npred);


    double *anggrid = (double*) malloc(sizeof(double) * Npred);
    vect_Linspace(anggrid, 0.0, 2 * PI, Npred);
    double *basgrid = (double*) malloc(sizeof(double) * Npred);
    vect_Linspace(basgrid, 4.3, 19.85, Npred);
    
    // intgrid is never used.
    //double *intgrid = (double*) malloc(sizeof(double) * Npred);
    
    
    

    int *inds = (int*) malloc(sizeof(int)* npts);
    vect_Sort(ang0Vect, inds, npts);
    double *xptemp = (double*) malloc(sizeof(double) * npts);
    for(i= 0; i<ndp; i++){
        xptemp[i] = ang0Vect[i];
    }
    qsort(xptemp, ndp, sizeof(double), compare);
    //for(i = 0; i<npts; i++) printf("%3.4f \n", xptemp[i]);
    double *yptemp = (double*) malloc(sizeof(double) * ndp);
    for(i = 0; i<ndp; i++){
        yptemp[i] = volumeVect[inds[i]];
    }
    
    
    double *xp = (double*) malloc(sizeof(double) * (ndp+2));
    double *yp = (double*) malloc(sizeof(double) * (ndp+2));
    xp[0] = xptemp[ndp-1] - (2 * PI);
    yp[0] = yptemp[ndp-1];
    for (i = 1; i < (ndp+1); i++){
        xp[i] = xptemp[i-1];
        yp[i] = yptemp[i-1];
    }
    xp[ndp+1] = xptemp[0] + (2*PI);
    yp[ndp+1] = yptemp[0];
    //vect_Print(yp, ndp+2, 4);
    
    
    int lenxp = ndp + 2;
    double vgridmax = vect_Max(vgrid, Npred);
    double vgridmin = vect_Min(vgrid, Npred); 
    //double **conhull = create2Darray(lenxp, Npred,vgridmax);
    //double **conhullm = create2Darray(lenxp, Npred, vgridmin);
    //matrix_Print4(conhull, lenxp, Npred, 3);
    //double **conhull, **conhullm;
    double *hullkeep = (double*) malloc(sizeof(double)*Npred);
    double *hullkeepm = (double*) malloc(sizeof(double) * Npred);
    for (i = 0; i< Npred; i++){
        hullkeep[i] = vgridmax;
        hullkeepm[i] = vgridmin;
    }
    double **conhull = create2Darray(lenxp, Npred,vgridmax);
	double **conhullm = create2Darray(lenxp, Npred, vgridmin);
	double *chultmp= (double*) malloc(sizeof(double)*lenxp);
	double *chulmtmp= (double*) malloc(sizeof(double)*lenxp);
	
	//vect_Print(hullkeepm, Npred, 4);
    int n, k, r;
    double m;


    for (n = 0; n<(lenxp -1); n++){
    	
		for(j = (n+1); j < lenxp; j++){
			 m = (yp[j] - yp[n])/(xp[j]-xp[n]);
			 //printf("(n, j):(%d, %d)  m = %3.4f \n",n, j,m);
			 
			for (k = 0; k < Npred; k++){
				//printf("(n, j):(%d, %d)   k = %d\n",n,j, k);
			 	if ((anggrid[k]<xp[n]) || (anggrid[k] > xp[j])){
			 		conhull[j][k] = vgridmax;
			 		//printf("option 1\n");
				}
				else{
			 		conhull[j][k] = m*(anggrid[k] - xp[n]) + yp[n];
			 		//printf("option 2\n");
			 	}
			 	if ((anggrid[k] <xp[n]) || (anggrid[k] > xp[j])){
			 		conhullm[j][k] = vgridmin;
			 		//printf("option 3\n");
			 	}
			 	else{
			 		conhullm[j][k] = m*(anggrid[k] - xp[n]) + yp[n];
			 		//printf("option 4\n");
			 	}
			}
		}
	
	
		for (k = 0; k<Npred; k++){
			for (r = 0; r< lenxp; r++){
				chultmp[r] = conhull[r][k];
				chulmtmp[r] = conhullm[r][k];
			}
			hullkeep[k] = min(hullkeep[k], vect_Min(chultmp, lenxp));
			hullkeepm[k] = max(hullkeepm[k],vect_Max(chulmtmp, lenxp));
		}
		//printf("hullkeepm[1] = %3.4f\n", hullkeepm[0]);
	}
	
	// No need for these matrices anymore.
	free(chultmp); free(chulmtmp);
	free2Darray(conhull, lenxp, Npred); free2Darray(conhullm, lenxp, Npred);
	
	//vect_Print(hullkeepm, Npred, 4);

   
    double thetas[] = {log(4), log(2), log(0.35)};
    int totsizeparams = (npts*6) + 10;
    double *wrparams = (double*) malloc(sizeof(double)*totsizeparams);
	// Wrapping up parameters into one array	
	wrapper(volumeVect, ang0Vect, basFricAngVect, intAngVect, 
           volBasVect, y, circterms, parest, npts, wrparams);
	
	/* !!! REMEMBER TO UNCOMMENT THIS !!!*/
	posteriorMode(thetas, wrparams);
	
	free(wrparams);
	//double output = wrprTimesLhood(thetas, wrparams); 
	//printf("out = %E\n", output);
	//thetas[0] = 2.2889; thetas[1] = 0.8176; thetas[2] = -1.1628;
    
    //double output = prTimesLhood (volumeVect, ang0Vect, basFricAngVect, 
                        //   intAngVect,volBasVect, y, thetas,
                        //    circterms, parest, npts);
	//printf("Output = %E\n", output);
   
 
    vect_Exp(thetas, thetas, 3);
    vect_Print(thetas, 3, 4);
    double theta1 = thetas[0];
    double theta2 = thetas[1];
    double theta3 = thetas[2];


    double alpha = 1.9;
    double alphat = 2.0;
   

	/* Correlation matrix */
    double **C = create2Darray(ndp, ndp, 0.0);
   	//matrix_PrintPart2(C, 0, 5, 0, 5);
    double *cn = (double*) malloc(sizeof(double) * 5);
	for (i= 0; i<5; i++){
		cn[i] = exp(-1 *pow(circterms[i], 2)/(4 * theta2))/
				sqrt(PI * theta2);
	}


    double c0 = 1.0/(2*sqrt(PI * theta2));
    int jj, kk;
    double *temp1 = (double*) malloc(sizeof(double) * 5);
    double corrangle;
    
    for (jj = 0; jj<ndp; jj++){
    	for(kk = 0; kk<ndp; kk++){
    		/* Calculate the correlation due to angle */
    		for(i = 0; i<5; i++){
    			temp1[i] = cn[i] * cos(circterms[i] * 
    									(ang0Vect[jj] - ang0Vect[kk]));
    		}
    	
    	corrangle = c0 + vect_Sum(temp1, 5);
    	C[jj][kk] = exp(-1*theta1 * 
    				pow_double(abs(volumeVect[jj] - volumeVect[kk]), alpha)) *
    				exp(-1*theta3 * 
    				pow_double(abs(basFricAngVect[jj] - 
    				basFricAngVect[kk]), alpha)) * corrangle;
    	}
    	
    }
	free(temp1); 
    
    double *Cinv = (double*) malloc(sizeof(double) * (ndp * ndp));
	int Cinv_index = 0;
	for(i = 0; i <ndp; i++){
		for (j = 0; j<ndp; j++){
			Cinv[Cinv_index] = C[j][i];
			Cinv_index++;
		}
	}
	matrix_Inverse(Cinv, ndp);
	free2Darray(C, ndp, ndp); // No longer need C as a 2D array. Cinv = C^-1
	
    double *XX = (double*) malloc(sizeof(double) * npts);
    double *YY_XX = (double*) malloc(sizeof(double) * ndp);
    

    for (i = 0; i<npts; i++){
    	XX[i] = beta_v + (m_v * volumeVect[i]) + (m_bf * basFricAngVect[i]) + 
    			(m_vbf * volBasVect[i]);
    	YY_XX[i] = y[i] - XX[i];
    	}
	
	double *logv = (double*) malloc(sizeof(double) * bofv_len);
	for(i = 0; i<bofv_len; i++){
		logv[i] = bofv_v[i] - 6;
	}


	//double *YY = (double*) malloc(sizeof(double) * npts);
	//for (i = 0; i<ndp; i++){
		//YY[i] = y[i];
    //}
    
    double *bweights = (double*) malloc(sizeof(double) * ndp);
    matrix_Mult('N', 'N', ndp, 1, ndp, 1.0, Cinv, ndp,YY_XX, 
                ndp, 0.0, bweights, ndp);
    
	//vect_Print(bweights, ndp, 4);
	/*bw2 = C\y = C\YY : YY = y*/
    double *bw2 = (double*) malloc(sizeof(double) * ndp);
    matrix_Mult('N', 'N', ndp, 1, ndp, 1.0, Cinv, ndp,y, 
                ndp, 0.0, bw2, ndp);

    
    double *bw3 = (double*) malloc(sizeof(double) * ndp);
	double *onestemp = (double*) malloc(sizeof(double) * ndp);
	for(i = 0; i<ndp; i++){
		onestemp[i] = 1.0;
	}
	matrix_Mult('N', 'N', ndp, 1, ndp, 1.0, Cinv, ndp,onestemp, 
                ndp, 0.0, bw3, ndp);
    free(onestemp);

    
    double d2 = 0;
    for(i = 0; i<(ndp * ndp); i++){
    	d2 += Cinv[i];
    }
	double n2 = 0;
	double n3 = 0.0;
	double *n3temp = (double *) malloc(sizeof(double) * ndp);
	matrix_Mult('N', 'N', ndp, 1, ndp, 1.0, Cinv, ndp,volumeVect, 
                ndp, 0.0, n3temp, ndp);
	for(i = 0; i<ndp; i++){
		n2 += bw2[i];
		n3 += n3temp[i];
	}
	free(n3temp);
	
	double VR = n3/d2;
	double *VV = (double *) malloc(sizeof(double) * ndp);
	for(i = 0; i<ndp; i++){
		VV[i] = volumeVect[i] - VR;
	}
	
	double n4 = 0.0;
	double d4 = 0.0;
	double *d4temp = (double *) malloc(sizeof(double) * ndp);
	matrix_Mult('N', 'N', ndp, 1, ndp, 1.0, Cinv, ndp,VV, 
                ndp, 0.0, d4temp, ndp);
	for(i = 0; i<ndp; i++){
		n4 += VV[i] * bw2[i];
		d4 += VV[i] * d4temp[i];
		
	}
	//printf("n4 = %3.5f\t d4= %3.4f\n",n4, d4);
	free(d4temp);
	double *bw4 = (double*) malloc(sizeof(double) * ndp);
	matrix_Mult('N', 'N', ndp, 1, ndp, 1.0, Cinv, ndp,VV, 
                ndp, 0.0, bw4, ndp);
    
	

	/* These curves are the volume-basal friction relationship */
    FILE *fLines;
    fLines = fopen("line_samples_Danilo.txt", "r");
    int nlines = linesInFile(fLines);
    printf("%d volume-basFric curves\n", nlines);
    double *linesamples1 = (double*) malloc(sizeof(double) * nlines);
    double *linesamples2 = (double*) malloc(sizeof(double) * nlines);			
    readLineSamples(fLines, nlines, linesamples1, linesamples2);
	fclose(fLines);
	//vect_Print(linesamples2, nlines, 4);
	
	

	
	double ***ypredfull = create3Darray2(Npred, Npred, Npred, 0.0);
	double ***ypredfullb = create3Darray2(Npred, Npred, Npred, 0.0);
	
    
    double *corrangles = (double*) malloc(sizeof(double) * Npred);
    double *temp = (double*) malloc(sizeof(double) * Npred);
    double *diffangle = (double*) malloc(sizeof(double) * Npred);
    
    int lind, nn, ii;
    
    double ypredtemp1, ypredtemp2;
    
   
    for (jj = 0; jj <ndp; jj++){
    	memset(corrangles, 0, sizeof(double) * Npred);
    	memset(temp, 0, sizeof(double) * Npred);
		
		/*for(i = 0; i<Npred; i++){
			//corrangles[i] = 0.0;
			//temp[i] = 0.0;
			diffangle[i] = anggrid[i] - ang0Vect[jj];
		}*/
		
		// diffangle = anggrid[i] - ang0Vect[jj] 
		for (lind = 0; lind < 5; lind++){
			for (i = 0; i<Npred; i++){
				temp[i] += cn[lind] * cos(circterms[lind] 
											* (anggrid[i] - ang0Vect[jj]));
			}
		}
		
		for (i = 0; i<Npred; i++){
			corrangles[i] = temp[i] + c0;
		}
		for (nn = 0; nn< Npred; nn++){
			for (kk=0; kk<Npred; kk++){
				
				ypredtemp1 = exp(-1*theta1 * pow_double(abs(vgrid[kk] -
								volumeVect[jj]), alpha));
				
				ypredtemp2 = exp(-1*theta3 *pow_double(abs(basgrid[nn] -
								basFricAngVect[jj]), alpha));
				/* Calculate the correlation due to angle ypred is mean 
				surface using MLE plugin*/
				for (i=0; i<Npred; i++){
					ypredfull[nn][kk][i] = ypredfull[nn][kk][i] + 
											(bweights[jj] * corrangles[i] *
											ypredtemp1 * ypredtemp2);
											
					ypredfullb[nn][kk][i] = ypredfullb[nn][kk][i] + 
											((bw2[jj] - ((n2*bw3[jj])/d2) - 
											((n4*bw4[jj])/d4))*corrangles[i]*
											ypredtemp1 * ypredtemp2);
					//printf("... %3.4f ", ypredfullb[nn][kk][i]);
					
				}
			
			}
		
		}
		//printf("ypredtemp2 = %3.4f\n", ypredtemp2);
	}
    

	/* Elts of ypredfull are in the order of 10^-12. use scientific notation */
    double ***ypredlbsave = create3Darray2(Npred, Npred, Npred, 0.0);
    matrix_Copy3Darray(ypredfullb, ypredlbsave, Npred, Npred, Npred);
	
	double ***ypredsave = create3Darray2(Npred, Npred, Npred, 0.0);
	matrix_Copy3Darray(ypredfull, ypredsave, Npred, Npred, Npred);
	
    
    double ***yreg = create3Darray2(Npred, Npred, Npred, 0.0);
    
    int mm;

    for(mm = 0; mm < Npred; mm++){
    	for(nn = 0; nn < Npred; nn++){
    		for (i = 0; i<Npred; i++){
    			yreg[i][mm][nn] = yreg[i][mm][nn] + beta_v + (m_v * vgrid[mm]);
    		}
    	}    	
    	for (i = 0; i<Npred; i++){
    		for(j= 0; j<Npred; j++){
    			ypredfullb[i][mm][j] = ypredfullb[i][mm][j] + 
    									((n2/d2)+((n4/d4)*(vgrid[mm]-VR)));
    		}
    	}
	}

	for (i = 0; i<Npred; i++){ // plane
    	for(j = 0; j<Npred; j++){ //row 
    		for(k = 0; k<Npred; k++){ //column
    			ypredfull[i][j][k] = ypredfull[i][j][k] + yreg[i][j][k];
    		}
    	}
    }


	// CREATE MESHGRID
	
/*	 
	// DO WE NEED TO READ IN DANILO_QS.TXT? We arent using it. NO..

	FILE *fDanilo_qs;
	fDanilo_qs = fopen("danilo_qs.txt", "r");
	int nlines_danilo_qs = linesInFile(fDanilo_qs);
	int ncols_danilo_qs = colsInFile(fDanilo_qs);
	double **danilo_qs = create2Darray(nlines_danilo_qs, ncols_danilo_qs, 0.0);
	readDaniloQs(fDanilo_qs, nlines_danilo_qs, danilo_qs);
    fclose(fDanilo_qs);
    
	double *q1 = (double*) malloc(sizeof(double) * nlines_danilo_qs);
	double *q99 = (double*) malloc(sizeof(double) * nlines_danilo_qs);
	
	for (i = 0; i<nlines_danilo_qs; i++){
		q1[i] = danilo_qs[i][5];
		q99[i] = danilo_qs[i][6];
	}
*/

	
	
  
    /* Read in quantiles. */
    FILE *fQuants;
    fQuants = fopen("quants.txt", "r");
    int nlines_quants = linesInFile(fQuants);
    //printf("nlinesquants =  %d\n", nlines_quants);
    double *vols = (double*) malloc(sizeof(double) * nlines_quants);
    double *q1 = (double*) malloc(sizeof(double) * nlines_quants);
    double *q5 = (double*) malloc(sizeof(double) * nlines_quants);
    double *q50 = (double*) malloc(sizeof(double) * nlines_quants);
    double *q95 = (double*) malloc(sizeof(double) * nlines_quants);
    double *q99 = (double*) malloc(sizeof(double) * nlines_quants);
	readQuants(fQuants, nlines_quants,vols,q1,q5,q50,q95,q99);
    fclose(fQuants);
    
    
    int hh;
    
    /* cutoff : ...*/
    double cutoff = 7.0; 
    int indrand;
    double *bofvtemplst = (double*) malloc(sizeof(double) * bofv_len);
    double *bofv = (double*) malloc(sizeof(double) * Npred); // Npred = 100

	double **ypred = create2Darray(Npred, Npred, 0.0);   		
   	double **ypredb = create2Darray(Npred, Npred, 0.0);
    double **ypredtemp = create2Darray(Npred, Npred, 0.0);
    double **ypredtempb = create2Darray(Npred, Npred, 0.0);
	double *zinterptemp = (double*) malloc(sizeof(double) * Npred);
	double *zinterptempb = (double*) malloc(sizeof(double) * Npred);
    double *alldist = (double*) malloc(sizeof(double) * Npred);
    
    double *vcontour = (double*) malloc(sizeof(double) * Npred);
    double *vcontourb = (double*) malloc(sizeof(double) * Npred);
    
    double yt;
    double dv= vgrid[Npred-1] - vgrid[0];
    int lb;
    double dist2h;
    double w;
    double thres = log(2);
    
    double *hz = (double*) malloc(sizeof(double) * Npred);
    double *signhz = (double*) malloc(sizeof(double) * (Npred-1));
    
    int zindsnbr; 
    int counter;
    int indp; int indm;
    
    double vcontourtemp;
    double maxypredcoltemp;
    
    int Nc = 50;
	double **vcontoursave = create2Darray(Nc, Npred, 0.0);
	double **vcontoursaveb = create2Darray(Nc, Npred, 0.0);
	double **bofvsave = create2Darray(Nc, Npred, 0.0);
	int *indrandsave = (int*) malloc(sizeof(int) * Nc);
    
    
    //vect_checkMonot(bofvtemplst, bofv_len);
    //vect_checkMonot(q1, nlines_quants);
    //vect_Print(q1, nlines_quants, 3);
    
    //vect_Print(logv,bofv_len, 4);

 //   #if 0
    for(hh = 0; hh< Nc; hh++){
    	
    	if(((hh%10) == 0) || (hh == (Nc-1))){
    		printf("hh = %d\n", hh);
    	}
    	
    	//indrand = hh;
    	indrand = (int)ceil((((double)rand())/(double)RAND_MAX) * 20000);
    	indrandsave[hh] = indrand;
   		for (i = 0; i<bofv_len; i++){
			bofvtemplst[i] = (180/PI) * atan(exp(linesamples1[indrand] + 
								(linesamples2[indrand] * logv[i])));

    	}
    	
    	if (hh > 4){
    		vect_Interpolation(bofv_v, bofvtemplst, vgrid, bofv, bofv_len, 
    							Npred);
    	}
    	else if (hh == 0){
    		vect_Interpolation(vols, q1, vgrid, bofv, bofv_len, Npred);
    	}
    	else if (hh == 1){
    		vect_Interpolation(vols, q5, vgrid, bofv, bofv_len, Npred);
    	}
    	
    	else if (hh == 2){
    		vect_Interpolation(vols, q50, vgrid, bofv, bofv_len, Npred);
    	}
    	else if (hh == 3){
    		vect_Interpolation(vols, q95, vgrid, bofv, bofv_len, Npred);
    		
		}
		else if (hh == 4){
			vect_Interpolation(vols, q99, vgrid, bofv, bofv_len, Npred);    	
    	}
    	
    	//vect_Print(bofv, Npred,4);
    	for (i = 0; i< Npred; i++){
    		if (bofv[i] < cutoff){
    			bofv[i] = 7+0.0001 * (double)rand() / (double)RAND_MAX; // rand;
    		}
    	}
    	
    	
    	for (i= 0; i< Npred; i++){
			bofvsave[hh][i] = bofv[i];
		}
		
		matrix_ZeroOut(ypred, Npred, Npred);
		matrix_ZeroOut(ypredb, Npred, Npred);

		for (k = 0; k < Npred; k++){
			
			for (i = 0; i<Npred; i++){
				for(j = 0; j<Npred; j++){
					ypredtemp[i][j] = ypredfull[i][j][k];
					ypredtempb[i][j] = ypredfullb[i][j][k];
				}
			}
			vect2DInterp(vgrid, basgrid, ypredtemp, vgrid, bofv,
						zinterptemp, Npred);
			vect2DInterp(vgrid, basgrid, ypredtempb, vgrid, bofv,
						zinterptempb, Npred);
			
			for (i = 0; i<Npred; i++){
				ypred[i][k] = zinterptemp[i];
				ypredb[i][k] = zinterptempb[i];
			}
		
		}

	
		/* ypredtemp = ypred. Not necessary since ypredtemp isnt used
		   again.*/

		yt = matrix_Min(ypred, Npred, Npred);
		/* dv is calculated outside of the for loop.
		 dist2hsave and wsave are neither filled up nor used. 
		 Will implement it later.
		*/ 
		//printf("beginning 2 for loops\n");
		for (k = 0; k<Npred; k++){
			
			lb = lookUpIndex(hullkeep[k], vgrid, Npred);
			//printf("lb = %d\n", lb);
			
			if (lb  != 0){
			//printf("lb = %d\n", lb);
			
			/* THE <= IS NECESSARY. IF LB is 1, has to do 2 iterations.*/
				for (j= 0; j<= lb; j++){
					
					for(i = 0; i<Npred; i++){
						alldist[i] = sqrt((pow(anggrid[k] - anggrid[i], 2)/
									pow(2*PI, 2)) + (pow(vgrid[j] - 
									hullkeep[i], 2)/pow(dv, 2)));
					}
				
					dist2h = vect_Min(alldist, Npred);
					//printf("dist2h = %3.4f\n", dist2h);
					
					w = exp((-1*pow(dist2h, 2)/2) * 14 * 14);
					//printf("dist2h=%3.4f\tw = %3.4f\n",dist2h, w);
					ypred[j][k] = w * ypred[j][k] + (1-w)*0*yreg[0][0][0];
					ypredb[j][k] = w * ypredb[j][k] + (1-w)*0*yreg[0][0][0];
				}
			}
		}
		
		
		/* ypred += 1e-9*randn(Npred, Npred);*/
		
		memset(vcontour, 0, sizeof(double) * Npred);
		memset(vcontourb, 0, sizeof(double) * Npred);
		
		for (jj = 0; jj < Npred; jj++){
			zindsnbr = 0;
			
			// create hz
			for (i = 0; i<Npred; i++){
				hz[i] = ypred[i][jj] - thres;
			}
			
			// create signhz & count how many are < 0
			for(i = 0; i<(Npred-1); i++){
				signhz[i] = hz[i] * hz[i+1];
				if (signhz[i] < 0){
					zindsnbr += 1;
				}
			}
			//printf("zindsnbr = %d\n", zindsnbr);
			
			// if zindsnbr is 0, fill vcontour with max(vgrid)
			if (zindsnbr == 0){
				vcontour[jj] = vect_Max(vgrid, Npred);
			}
			
			// if zindsnbr != 0, then create zindslist and fill vcontour 
			// accordingly
			else{
				int *zindsList = (int*) malloc(sizeof(int) * zindsnbr);
				counter = 0;
				// create the list of the indices of zinds
				for(i = 0; i<(Npred-1); i++){
					if (signhz[i] < 0){
						zindsList[counter] = i;
						counter++;
					}	
				}
				
				// If only one index
				if (zindsnbr == 1){
					indp = min(zindsList[0]+ 1, Npred-1);
					indm = max(zindsList[0], 0);
					//printf("indp: %d \t indm = %d\n", indp, indm);
					double *ypredcoltemp = (double*) malloc(sizeof(double) * 
											(indp - indm +1));
					double *vgridtemp = (double*) malloc(sizeof(double) * 
											(indp - indm +1));
					
					// counter = length of ypredcoltemp and vgridtemp
					counter = 0; 
					
					for(i = indm; i<= indp; i++){
						ypredcoltemp[counter] = ypred[i][jj];
						vgridtemp[counter] = vgrid[i];
						counter++; //at the end of forloop == indp - indm + 1
					}
					
					
					
					vect_Interpolation(ypredcoltemp, vgridtemp, &thres, 
										&vcontourtemp, counter, 1);

										 
					vcontour[jj] = vcontourtemp;
					
					free(ypredcoltemp);
					free(vgridtemp);
				}
				
				// If there are more than 1 index
				
				else{
					indp = min(zindsList[zindsnbr-1]+ 1, Npred-1);
					indm = max(zindsList[0], 0);
					
					/* variable-sized vectors. That's why they are declared
					   and freed inside these loops*/
					double *ypredcoltemp = (double*) malloc(sizeof(double) * 
											(indp - indm +1));
					double *vgridtemp = (double*) malloc(sizeof(double) * 
											(indp - indm +1));
					//counter = length of ypredcoltemp and vgridtemp
					counter = 0; 
					
					for(i = indm; i<= indp; i++){
						ypredcoltemp[counter] = ypred[i][jj];
						vgridtemp[counter] = vgrid[i];
						counter++; //at the end of forloop == indp - indm + 1
					}
				
					maxypredcoltemp = vect_Max(ypredcoltemp, counter);
					
					if (maxypredcoltemp > thres){
						vect_Interpolation(ypredcoltemp, vgridtemp, &thres, 
										&vcontourtemp, counter, 1);
						vcontour[jj] = vcontourtemp;
					}
					else{	
						vcontour[jj] = vect_Max(vgrid, Npred);
					}

					free(ypredcoltemp);
					free(vgridtemp);
				}
				free(zindsList);
			}
		}
	
		/* ===================================================================
		 =====================================================================*/
		
		
		for (jj = 0; jj < Npred; jj++){
			zindsnbr = 0;
			
			// create hz
			for (i = 0; i<Npred; i++){
				hz[i] = ypredb[i][jj] - thres;
			}
			
			// create signhz & count how many are < 0
			for(i = 0; i<(Npred-1); i++){
				signhz[i] = hz[i] * hz[i+1];
				if (signhz[i] < 0){
					zindsnbr += 1;
				}
			}

			
			// if zindsnbr is 0, fill vcontour with max(vgrid)
			if (zindsnbr == 0){
				vcontourb[jj] = vect_Max(vgrid, Npred);
			}
			
			// if zindsnbr != 0, then create zindslist and fill vcontour 
			// accordingly
			else{
				int *zindsList = (int*) malloc(sizeof(int) * zindsnbr);
				counter = 0;
				// create the list of the indices of zinds
				for(i = 0; i<(Npred-1); i++){
					if (signhz[i] < 0){
						zindsList[counter] = i;
						counter++;
					}	
				}
				
				// If only one index
				if (zindsnbr == 1){
					indp = min(zindsList[0]+ 1, Npred-1);
					indm = max(zindsList[0], 0);
					//printf("indp: %d \t indm = %d\n", indp, indm);
					double *ypredcoltemp = (double*) malloc(sizeof(double) * 
											(indp - indm +1));
					double *vgridtemp = (double*) malloc(sizeof(double) * 
											(indp - indm +1));
					
					// counter = length of ypredcoltemp and vgridtemp
					counter = 0; 
					
					for(i = indm; i<= indp; i++){
						ypredcoltemp[counter] = ypredb[i][jj];
						vgridtemp[counter] = vgrid[i];
						counter++; //at the end of forloop == indp - indm + 1
					}
					
					
					
					vect_Interpolation(ypredcoltemp, vgridtemp, &thres, 
										&vcontourtemp, counter, 1);
															 
					vcontourb[jj] = vcontourtemp;
					
					free(ypredcoltemp);
					free(vgridtemp);
				}
				
				// If there are more than 1 index
				
				else{
					indp = min(zindsList[zindsnbr-1]+ 1, Npred-1);
					indm = max(zindsList[0], 0);
					double *ypredcoltemp = (double*) malloc(sizeof(double) * 
											(indp - indm +1));
					double *vgridtemp = (double*) malloc(sizeof(double) * 
											(indp - indm +1));
					//counter = length of ypredcoltemp and vgridtemp
					counter = 0; 
					
					for(i = indm; i<= indp; i++){
						ypredcoltemp[counter] = ypredb[i][jj];
						vgridtemp[counter] = vgrid[i];
						counter++; //at the end of forloop == indp - indm + 1
					}
				
					maxypredcoltemp = vect_Max(ypredcoltemp, counter);
					
					if (maxypredcoltemp > thres){
						vect_Interpolation(ypredcoltemp, vgridtemp, &thres, 
										&vcontourtemp, counter, 1);
						vcontourb[jj] = vcontourtemp;
					}
					else{	
						vcontourb[jj] = vect_Max(vgrid, Npred);
					}
											
					
					
					
					free(ypredcoltemp);
					free(vgridtemp);
				}
				free(zindsList);
			}
		}

		/* Save contour lines*/
		for(i = 0; i<Npred; i++){
			vcontoursave[hh][i] = vcontour[i];
			vcontoursaveb[hh][i] = vcontourb[i];
		}	
    }
    
    

    clock_t toc = clock();
	printf("End of Run: %f secs\n", (double)(toc - tic) / CLOCKS_PER_SEC);
    
    
    char outdirname[MAXNAME];
    strcpy(outdirname, "./contours");
    FILE *fvcontoursave;
    char outfname[MAXNAME];
    snprintf(outfname, MAXNAME, "%s/contours40_%06d.txt", outdirname,
                idfile);
    puts(outfname);
    fvcontoursave = fopen(outfname, "w");

    matrix_PrintToFileDOUBLE(vcontoursaveb, fvcontoursave, Nc, Npred, 4);
    fclose(fvcontoursave);
    
    
    FILE *fvgrid;
    char vgriddir[MAXNAME];
    strcpy(vgriddir, "./vgrids");
    snprintf(outfname, MAXNAME, "%s/vgrid_%06d.txt",vgriddir,
                idfile);
    fvgrid = fopen(outfname, "w");
    vect_PrintToFileDOUBLE(vgrid, fvgrid, Npred, 4);
    fclose(fvgrid);
    /*FILE *fanggrid;
    fanggrid = fopen("anggrid.txt", "w");
    vect_PrintToFileDOUBLE(anggrid, fanggrid, Npred, 4);
    fclose(fanggrid);
    */


    FILE *fypredb;
    char ypredbdir[MAXNAME];
    strcpy(ypredbdir, "./ypredbs");
    snprintf(outfname, MAXNAME, "%s/ypredb_%06d.txt", ypredbdir, idfile);
    fypredb = fopen(outfname, "w");
    matrix_PrintToFileDOUBLE(ypredb, fypredb, Npred, Npred, 4); 
    fclose(fypredb);
    //matrix_PrintPart3(ypredb, 0, 4, 0, 4);
    
	
	//#endif


    return 0;


}

