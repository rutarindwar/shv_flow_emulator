#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <f2c.h>
#include <blaswrap.h>
#include <clapack.h>


#define MAXCOLS 100
#define MAXROWS 500


void vect_Exp(double *lst, double *result, int len){

    int i;
    for (i = 0; i < len; i++){
        result[i] = exp(lst[i]);
    }
}

void vect_checkMonot(double *vect, int len){
	int i;
	int flag = 1;
	for(i = 1; i<len; i++){
		if (! (vect[i-1] < vect[i])){
			flag = 0;
			break;
		}
	}
	if (flag){
		printf("Monotonic\n");
	}
	else{
		printf("NOT MONOTONIC (%d)!!!!\n", i-1);
	}
	
}	


void vect_Print(double *lst, int len, int precision){
     /* Prints out an array of npts doubles*/

     int j;
     printf("[");
     for(j = 0; j<len ; j++){
           printf(" %4.*f", precision, lst[j]);
     }
     printf("]\n");
}

void vect_PrintINT(int *lst, int len){
     /* Prints out an array of npts doubles*/

     int j;
     printf("[");
     for(j = 0; j<len ; j++){
           printf(" %d", lst[j]);
     }
     printf("]\n");
}


void vect_PrintPart(double *lst, int a, int b, int precision){
    int i;

    for (i = a; i<= b; i++){
        printf("%3.*f ",precision, lst[i]);
    }
    printf("\n");
}



void vect_MultN(double *lst, double *result, double k, int len){

    int i;
    for (i = 0; i<len; i++){
        result[i] = lst[i] * k;
    }
}


void matrix_Print1(double* mat, int nrows, int ncols){

    int i;
    int totsize = nrows * ncols;
    printf("number of entries: %i\n", totsize);
    for (i= 0; i<(totsize); i++){
        if (((i+ 1)% ncols) == 0) {
            printf("%3.5f\n", mat[i]);
        }
        else{
            printf("%3.5f ", mat[i]);
        }
    }
}



void matrix_Print4(double **mat, int nrows, int ncols, int precision){
/*Print a 2D array created by using malloc twice. 
 * Need to give it a pointer to a pointer.
 use matrix created with create2Darray*/
    int i, j;
    for(i= 0; i<nrows; i++){
        for(j = 0; j <ncols; j++){
            printf(" %3.*f", precision, mat[i][j]);
        }
        printf("\n");
    }
} 


double ***create3Darray(int nrows, int ncols, int nplanes, double val){

	double ***array3D = (double***) malloc(sizeof(double**) * nrows);
	
	int i, j,k;
	i = 0; j = 0; k = 0;
	for(i = 0; i<nrows; i++){
		array3D[i] = (double**) malloc(sizeof(double*) * ncols);
		//printf("created %dth column\n", i+1);
		for (j = 0; j < ncols; j++){
			array3D[i][j] = (double*) malloc(sizeof(double) * nplanes);
			//printf("created %dth plane\n", j+1);
		}
	}

	return array3D;
		
}



double **create2Darray(int nrows, int ncols, double val){
/* Create a 2d array using malloc. Can be indexed with [i][j].
 * To be used along with matrix_Print4 */
    double **array2D = (double**) malloc(sizeof(double*) * nrows);
    //printf("created array2d of nrows: %d  ncols: %d\n", nrows, ncols);
    int i; int j;
    for(i= 0; i< nrows; i++){
        array2D[i] = (double*) malloc(sizeof(double)*ncols);
        //printf("created %dth column\n", i+1);
        for (j = 0; j< ncols; j++){
            array2D[i][j] = val;
            //printf("%3.4f ", array2D[i][j]);
        }
        //printf("\n");
    }
    return array2D;
}


void free2Darray(double **mat, int nrows, int ncols){
    int i;
    for(i = 0; i<nrows; i++){
        free(mat[i]);
    }
    free(mat);
}


double ***create3Darray2(int nplanes, int nrows, int ncols, double val){

	double ***array3D = (double***) malloc(sizeof(double**) * nplanes);
	
	int i;
	for(i = 0; i <nplanes; i++){
		array3D[i] = create2Darray(nrows, ncols, val);
	}
	
	return array3D;	
}





void matrix_Print2 (double* mat,int nrows, int ncols, int precision){
    /*myPrintArray(&matrix[0][0], r, c);*/

    int totsize = nrows * ncols;

    int i;
    for (i = 0; i< totsize; i++){
            if (((i+1)%ncols) == 0){
                printf("%3.*f\n", precision,*(mat+i));
            }
            else{
                printf("%3.*f ", precision, *(mat+i));
            }
    }
    printf("\n");
}



void matrix_MultN(double* mat, double k, int nrows, int ncols){

    int totsize = nrows*ncols;
    int i;
    for (i = 0; i< totsize; i++){
            *(mat + i) = *(mat+i) * k;
    }
}

void vect_Power(double *lst, double *result, int k, int len){

    int i;
    for(i = 0; i<len; i++){
       result[i] = pow(lst[i], (k*1.0));
    }
}


double vect_Sum(double *lst, int len){

    double sum = 0.0;
    int i;
    //vect_Print(lst, 5);
    for (i = 0; i <len; i++){
        sum = sum + lst[i];
    }
    //printf("Inside:    %3.4f \n", sum);
    return (sum);
}


void matrix_PrintPart(double* mat, int nrows, int ncols, int ri,
                      int rj, int ci, int cj, int precision){

    int i,j;
    double* b0 = mat + ci + (ncols*ri);
    double* b = b0;
    //printf("beginning: %3.4f\n\n", *b0);

    for (i = 0; i <= (rj-ri); i++){
        for (j = 0; j<= (cj-ci); j++){

            printf("%3.*f ",precision, *(b + j));
        }
        printf("\n");
        b = b0 + (ncols*(i+1));
    }
}

void matrix_PrintPart2(double **mat, int c0, int c1, int r0, int r1,
						int precision){
	/* Print a part of a 2d matrix. This matrix has to be the kind that is generated by create2Darray ( a pointer to a pointer of a double)
	!! PRINTS IN FLOATING POINT NOTATION !!*/
	
	int i,j;
	for(i = c0; i<= c1; i++){
		for(j= r0; j<= r1; j++){
			printf(" %3.*f", precision, mat[i][j]);
		}
		printf("\n");
	}
}

void matrix_PrintPart3(double **mat, int c0, int c1, int r0, int r1){
	/* Print a part of a 2d matrix. This matrix has to be the kind that is generated by create2Darray ( a pointer to a pointer of a double)
	!! PRINTS IN SCIENTIFIC NOTATION !!*/
	
	int i,j;
	for(i = c0; i<= c1; i++){
		for(j= r0; j<= r1; j++){
			printf(" %E", mat[i][j]);
		}
		printf("\n");
	}
}




void vect_AddN(double lst[], double res[],  int len, double k){

    int i;
    for (i = 0; i<len; i++){
        res[i] = lst[i] + k;
    }
}

        

void matrix_Print3 (double* mat,int nrows, int ncols, int precision){
    /*myPrintArray(&matrix[0][0], r, c);*/

    int totsize = nrows * ncols;
    int i, j;
    double *ind; ind = mat;
    for (i= 0; i<nrows; i++){
        for(j = 0; j<ncols; j++){
            printf("%3.*f ", precision, *(ind +(j* nrows)));
        }
        printf("\b\n");
        
        ind++;
    }
    printf("**************************************************\n");
}



void matrix_Inverse(double* mat, int size){

    long int m = size;
    long int n = size;
    long int lda = size;
    
    long int *ipiv;
    ipiv = (long int*) malloc(sizeof(long int) * size);
    long int info;

    dgetrf_(&m, &n,mat, &lda, ipiv, &info);
    
    long int lwork = size*size;
    double work[lwork];

    dgetri_(&n, mat, &lda, ipiv, work, &lwork, &info);

}


void matrix_sqrt(double* mat, int size, double* sqmat){

    char JOBZ = 'V'; char UPLO = 'U';
    long int n = size;
    long int lda = size;
    double w[size];
    long int lwork = (3*n) - 1;
    double work[lwork];
    
    long int INFO;

    dsyev_(&JOBZ, &UPLO, &n, mat,&lda, w, work, &lwork, &INFO);
    
   // vect_Print(w, 3,4);
    int totsize = size*size;
    double* diageig;
    diageig = (double*)  malloc(sizeof(double) * totsize);
    
    int i; double eigtemp;
    int eigind = 0;
    for (i = 0; i < totsize; i++){
        if (i %(size+1) == 0){
            eigtemp = w[eigind];
            diageig[i] = sqrt(eigtemp);
            eigind++;
        }
        else{
            diageig[i] = 0.0;
        }
    }

    char TRANSA = 'N';
    char TRANSB = 'N';
    long int nm = size, nn = size, nk = size;
    double alpha = 1.0;
    long int nlda = size;
    long int ldb = size;
    double nbeta = 0.0;
    double* tempC;
    tempC = (double*)  malloc(sizeof(double) * totsize);

    long int nldc = size;
    dgemm_(&TRANSA, &TRANSB, &nm, &nn, &nk, &alpha, mat, &nlda, 
           diageig, &ldb, &nbeta, tempC, &nldc);

    char TRANSC = 'T';
    dgemm_(&TRANSA, &TRANSC, &nm, &nn, &nk, &alpha, tempC, &nlda, 
            mat, &ldb, &nbeta, sqmat, &nldc);

    //vect_Print(w, size, 4);
   // matrix_Print3(sqmat, size, size, 4);

} 


void matrix_Mult(char tra, char trb, int m, int n, int k,
                 double al, double *mata, int lda, 
                 double *matb, int ldb, double bet, 
                 double *matc, int ldc){

    char TRANSA = tra, TRANSB = trb;
    long int nm = m, nn = n, nk = k, nlda = lda, nldb = ldb,
             nldc = ldc; 
    double nalpha = al, nbeta = bet;    

    dgemm_(&TRANSA, &TRANSB, &nm, &nn, &nk, &nalpha, mata,
            &nlda, matb, &nldb, &nbeta, matc, &nldc);
}  


void matrix_det(double *mat, int nrows, int ncols, double *det){

    long int m = nrows, n = ncols, lda = nrows, INFO;
    long int *ipiv = (long int*) malloc(sizeof(long int) *
                      min(nrows, ncols));
    
    dgetrf_(&m, &n, mat,&lda, ipiv, &INFO);

    int i;
    (*det) = 1;
    int size = ncols;
    int totsize = size*size;
    int count = 0;
    
    for(i = 0; i<totsize; i+=(size+1)){
                 
        if (ipiv[count] == (count+1)){
            (*det) = (*det) * (abs(mat[i]));
        }
        else{
            (*det) =(-1) *  (*det) * (abs(mat[i]));
        } 
        count++;
    }
    //printf("det = %E\n", (*det));
}


void matrix_PrintEnd(double *mat, double size, int len){

    int i;
    printf("...\n");
    for(i = (size-len); i<size; i++){
        printf("%3.4f\n", mat[i]);
    }
    printf("=========\n");
}


void matrix_Trace(double* mat, double *trace, int size){

    int totsize = size * size;
    int i;
    (*trace) = 0;
    //printf("size: %d\n", size);
    for(i = 0; i<totsize; i+=(size+1)){
        (*trace) = (*trace) + (mat[i]);
    }
}
 


void vect_Prod(double *lst, double *prod, int len){

    (*prod) = 1;
    int i;
    for (i = 0; i<len; i++){
        (*prod) = (*prod) * lst[i];
    }

}

/*
void matrix_solveAxb(double *a, double *b,  int size){

    long int nm = size;
    long int nn = size;
    long int lda = size;
    long int *ipiv;
    ipiv = (long int*) malloc(sizeof(long int) * size);
    long int info;




    dgetrf_(&nm, &nn,a, &lda, ipiv, &info);
    long int nrhs = 1; long int ldb = size;
    matrix_Print3(a, 2, 2) 
    dgesv_(&nm, &nrhs, a, &nn, ipiv, b, &ldb, &info);
    printf("Inside solveaxb, solution: \n"); 
    vect_Print(b, size, 4);
}

*/

void vect_LinFit (double *x, double *y, int len, double *coef){
    
    double *atemp = (double*) malloc(sizeof(double) * (len*2));
    int i;
    for (i = 0; i<len; i++){
        atemp[i] = 1.0;
    }
    for (i = len; i < (len*2); i++){
        atemp[i] = x[(i-len)];
    }
   // matrix_Print3(atemp, len, 2, 4);
    double *ata = (double*) malloc(sizeof(double)*4);
   
    long int m = 2, n = 2, k = len;
    double alpha = 1.0;
    long int lda = len, ldb = len;
    double beta = 0.0;
    long int ldc = 2;
    char TRANSA = 'T'; char TRANSB = 'N';
    dgemm_(&TRANSA, &TRANSB, &m, &n, &k, &alpha, atemp, &lda, atemp,
           &ldb, &beta, ata, &ldc); 
    
    long int nn = 1;
    double *atb = (double*) malloc(sizeof(double) * 2);
    dgemm_(&TRANSA, &TRANSB, &m, &nn, &k, &alpha, atemp, &lda, y, 
           &ldb, &beta, atb, &ldc);
    
   // matrix_Print3(ata, 2, 2, 4);
   // matrix_Print3(atb, 2, 1, 4);
    

    //double *sol = (double*) malloc(sizeof(double) * len);
    matrix_Inverse(ata, 2);
    char TRANSA2 = 'N';
    dgemm_(&TRANSA2, &TRANSB, &m, &nn,&n, &alpha, ata, &m, 
            atb, &n, &beta,coef, &m); 

    //vect_Print(sol, 2, 4);
}



void vect_Linspace(double *lst, double a, double b, int len){
    
    double step = (b-a)/(len-1);
    int i;
    double temp = a;
    for (i = 0; i<len; i++){
        lst[i] = temp;
        temp += step;
    }
    //vect_Print(lst, len,4);

}

double vect_Max(double *lst, int len){

    double maxm = lst[0];
    int i;
    for (i = 1; i<len; i++){
        if (lst[i] > maxm){
            maxm = lst[i];
        }
    }
    return maxm;
}

double vect_Min(double *lst, int len){
    double minm = lst[0];
    int i;
    for (i = 1; i<len; i++){
        if (lst[i] < minm){
            minm = lst[i];
        }
    }
    return minm;
}

struct pair{
    double a;
    int b;
};


int compare (const void *const first, const void *const second){
    
    if (((const struct pair *)first)->a >
        ((const struct pair *)second)->a){
        return 1;
    }
    else if (((const struct pair *)first)->a <
             ((const struct pair *)second)->a){
        return -1;
    }
    else
        return 0;
}





/*


int compare (const void * a, const void * b)
{
    if ( *(double*)a <  *(double*)b ) return -1;
    if ( *(double*)a == *(double*)b ) return 0;
    if ( *(double*)a >  *(double*)b ) return 1;    
}


*/



void vect_Sort(double *lst, int *inds, int len){

/* Sort a vector and return the indices of the sorted elements */

    double *lst_copy = (double*) malloc(sizeof(double) * len);
    int i;
    for (i = 0; i<len; i++){
        lst_copy[i] = lst[i];
    }
   
    

    qsort(lst_copy, len, sizeof(double), compare);
    
    int j;

    for (i = 0; i<len; i++){
        for (j = 0; j<len; j++){
            
            if (lst_copy[i] == lst[j]){
                inds[i] = j;
                break;
            }
        }
    }
    
//    for (i = 0; i<len; i++){
}
                
                        
void matrix_Copy3Darray(double ***source, double ***dest, int nr, int nc, 
						int np){
						
	
	int i, j,k;
	for (i= 0; i<np; i++){
		for(j= 0; j<nr; j++){
			for(k = 0; k < nc; k++){
				dest[i][j][k] = source[i][j][k];
			}
		}
	}
	//printf("Done copying 3D array\n");
}






void matrix_Transpose(double **mat, int nrows, int ncols){
	
	double temp;
	int i,j;
	for(i = 0; i<nrows; i++){
		for(j= 0; j<i; j++){
				printf("swapping .. \n");
				temp = mat[i][j];
				mat[i][j] = mat[j][i];
				mat[j][i] = temp;
		}
	}

}

void matrix_ZeroOut(double **mat, int nrows, int ncols){
	
	int i,j;
	
	for (i = 0; i<nrows; i++){
		memset(mat[i], 0, sizeof(double) * ncols);
	}
	
}



int lookUpIndex(double x, double *foo, int len){
	/* Goes with vect2DInterp. Search a grid and return index values to be
	used in vect2DInterp. Keep track of edge values.*/
	int ind, val;
    for (ind = 0; ind < len; ind++){
		if ((x == foo[ind]) && (ind == 0)){
			val = ind;
			//printf("found it & first elt\n");
			break;
		}
		if ((x == foo[ind]) && (ind == (len-1))){
			val = ind -1;
			//printf("last elt\n");
			break;
		}
		
		if ((x == foo[ind]) && (ind != 0) && (ind != (len-1))){
			val = ind;
			//printf("found it. neither first nor last elt\n");
			break;
		}
		
		else{
			if (x < foo[ind]){
				val = ind - 1;
				//printf("inside elt\n");
				break;
			}
		}
	}
    //printf("index = %d\n", val);
	return val;
}
			


void vect2DInterp(double *xvect, double *yvect, double **zgrid, double *xi,
					double *yi, double *zi, int len){
	/* Return a vector of interpolated values at the values (x,y) 
	   in xi and yi.
	   Modify this to take in the lengths of BOTH xvect and xi*/				
					

	int i;
	int xind, yind; 
	double z00, z01, z10, z11;
	double xm, ym;
	double zc1, zc2, zcm;
	for(i = 0; i < len; i++){
		
		xind = lookUpIndex(xi[i], xvect, len);
	//	printf("%3.4f\t < %3.4f\t < %3.4f\n", x[xind], xi[i], x[xind+1]);
	
		yind = lookUpIndex(yi[i], yvect, len);
		
		z00 = zgrid[yind][xind];
		z01 = zgrid[yind][xind +1];
		z10 = zgrid[yind+1][xind];
		z11 = zgrid[yind+1][xind+1];
		
		xm = xi[i];
		ym = yi[i];
		
		zc1 = z10 + (z11 - z10) * ((xm - xvect[xind])/(xvect[xind+1] - 
														xvect[xind]));
		zc2 = z00 + (z01 - z00) * ((xm - xvect[xind])/(xvect[xind+1] -
														xvect[xind]));
		
		zcm = zc2 + (zc1 - zc2) * ((ym - yvect[yind])/(yvect[yind+1] -
														yvect[yind]));
		
		zi[i] = zcm;
	
	}
	
	

}



void vect2DInterp2(double *xvect, double *yvect, double **zgrid,
                    double *xi, double *yi, double **zi, 
                    int xvlen, int yvlen, int xilen, int yilen){
	/* Return a grid of interpolated values at the values (x,y) 
	   in xi and yi. */				
					
	printf("Interpolating over grid ... ");
	int i = 0;
	int j = 0;
	int xind, yind; 
	double z00, z01, z10, z11;
	double xm, ym;
	double zc1, zc2, zcm;
	
	// Get boundaries to check for values outside the zgrid
	double xvectmin = vect_Min(xvect, xvlen);
	double xvectmax = vect_Max(xvect, xvlen);
	double yvectmin = vect_Min(yvect, yvlen);
	double yvectmax = vect_Max(yvect, yvlen);
    
    //printf("%e \t %e \n %e \t %e\n", xvectmin, xvectmax,
            //yvectmin, yvectmax);
    
	//printf("%d --- %d\n", xilen, yilen);
	//printf("%4.5f --- %4.5f\n", xvectmin, xvectmax);

    

    int count = 0;
    
	/*Note how I'm indexing zgrid, zgrid[j][i], and NOT [i][j]*/
	
	for(i = 0; i < yilen; i++){
        for(j = 0; j<xilen; j++){
            if ((yi[i] < yvectmin) || (yi[i] > yvectmax) ||
                (xi[j] < xvectmin) || (xi[j] > xvectmax)){
                    
					/*RANDOM PLACE HOLDER, used in place of 0, or NaN*/
					zi[j][i]= 0.0;
					
					
            //        printf("%e ", yi[i]);
              //      printf("%e", xi[j]);
                    //printf("(%d, %d) \t %d %d %d %d\n",i,j,
                           // (yi[i] < yvectmin), (yi[i] > yvectmax),
                           // (xi[j] < xvectmax), (xi[j] > xvectmax));                // printf("(%d, %d)\n", i, j);
                    count++;
            }
            else{
                
                //printf("----- (%d, %d)\n", i,j);
                xind = lookUpIndex(xi[j], xvect, xvlen);
                yind = lookUpIndex(yi[i], yvect, yvlen);

                //printf("xind = %d \t yind = %d\n", xind, yind); 
                
    			z00 = zgrid[yind][xind];
				//printf("z00 = %f \t ", z00);
				
                z01 = zgrid[yind][xind +1];
				//printf("z01 = %f \t", z01);
				
				
				z10 = zgrid[yind+1][xind];
				//printf("z10 = %f\n", z10);
				
				z11 = zgrid[yind+1][xind+1];
                
				//printf("%e \t %e \t %e \t%e\n",z00, z01, z10, z11); 
				
                xm = xi[j];
				ym = yi[i];
		
				zc1 = z10 + (z11 - z10) * ((xm - xvect[xind])/
									(xvect[xind+1] - xvect[xind]));
				zc2 = z00 + (z01 - z00) * ((xm - xvect[xind])/
									(xvect[xind+1] - xvect[xind]));
				zcm = zc2 + (zc1 - zc2) * ((ym - yvect[yind])/
									(yvect[yind+1] - yvect[yind]));
				zi[j][i] = zcm;
                
            }
        }
    }
    printf("Done\n");
    printf("NaN = %d\n", count);
    
}
















double matrix_Min(double **mat, int nrows, int ncols){
	
	int i, j;
	
	double mmin = mat[0][0];
	
	for (i = 0; i< nrows; i++){
		for(j = 0; j<ncols; j++){
			
			if (mat[i][j] < mmin){
				mmin = mat[i][j];
			}
		}
	}
	return mmin;
}




void matrix_PrintToFileDOUBLE(double **mat, FILE *fout, int nrows, 
								int ncols, int precision){
	
	int i;
	int j;
	//printf("Writing to file ...\n");
	for (i= 0; i< nrows; i++){
		for(j = 0; j < (ncols-1); j++){
			fprintf(fout, "%3.*f\t", precision, mat[i][j]);
		}
		fprintf(fout, "%3.*f\n", precision, mat[i][j]);
	}
	//printf("Done\n");
}


void vect_PrintToFileDOUBLE(double *vect, FILE *fout, int nrows, 
								int precision){
	int i;	
	for (i = 0; i<nrows; i++){
		fprintf(fout, "%3.*f\n", precision, vect[i]);
	}
}


int findIndex(double x, double *xlist, int len, int *xindex){
	
	int ind;
	int val;
	int found = 0;
	
    for (ind = 0; ind < len; ind++){
		if ((x == xlist[ind]) && (ind == 0)){
			found = 1;
			val = ind;
			//printf("found it & first elt\n");
			break;
		}
		if ((x == xlist[ind]) && (ind == (len-1))){
			found = 1;
			val = ind;
			//printf("last elt\n");
			break;
		}
		
		if ((x == xlist[ind]) && (ind != 0) && (ind != (len-1))){
			found = 1;
			val = ind;
			//printf("found it. neither first nor last elt\n");
			break;
		}
		
		else{
			if (x < xlist[ind]){
				found = 0;
				val = ind - 1;
				//printf("inside elt\n");
				break;
			}
		}
	}
    //printf("index = %d\n", val);
	(*xindex) = val;
	return found;
	
}


void vect_Interpolation(double *xvect, double *yvect, double *xi, double *yi,
							int xvect_len, int xi_len){

	//printf(" In interpolation2 \n");
	int i;
	int indtemp;
	int foundflag;
	double x0, x1, y0, y1;
	double xm, ym;
	for (i = 0; i<xi_len; i++){
		foundflag = findIndex(xi[i], xvect, xvect_len, &indtemp);
		
		if (foundflag){
			yi[i] = yvect[indtemp];
		}
		else{
			x0 = xvect[indtemp];
			y0 = yvect[indtemp];
			x1 = xvect[indtemp+1];
			y1 = yvect[indtemp+1];
			
			xm = xi[i];
			ym = y0 + (y1-y0) * ((xm-x0)/(x1-x0));
			yi[i] = ym;
		}
	}
	//printf("Out of interpolation2\n");
}
					




double matrix_Max(double **mat, int nrows, int ncols){
	/*Finds the maximum value in a 2d matrix*/
	double maxm = mat[0][0];
	int i,j;
	
	for(i =0 ; i<nrows; i++){
		for(j = 0; j<ncols; j++){
			if (mat[i][j] < maxm){
				maxm = mat[i][j];
			}
		}
	}
	
	return maxm;
}


double matrix_Sum(double **mat, int nrows, int ncols){
	/* Compute the total sum of a 2d matrix*/
	double msum = 0;	
	int i,j;
	for(i = 0; i<nrows; i++){
		for(j = 0; j<ncols; j++){
			msum += mat[i][j];
		}
	}
	return msum;
}
































