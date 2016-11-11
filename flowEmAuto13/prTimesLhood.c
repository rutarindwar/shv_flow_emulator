#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#include <gsl/gsl_poly.h>
//#include <gsl/gsl_complex.h>
//#include <gsl/gsl_spline.h>
//#include <gsl/gsl_math.h>

#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>




#include "matrixOperations.h"
#include "f2c.h"



double abs_double(double x){
    /* Return the Absolute value of a double number*/
    if (x>=0) return x;
    else return (-1*x);
}



double pow_double(double x, double y){
    /* Return x^y where x and y are doubles. This is useful in
     case x is negative and y is a double. (This would cause
     pow to take a n-root of a negative number,.. which returns #INF
    */

    int xsign;

    if (x>= 0){
        xsign = 1;
    }
    else{
        xsign = -1;
    }

    double temppow = pow(abs_double(x), y);

    if (fmod(y, 2.0)== 0){
        return temppow;
    }
    else{
        return xsign*temppow;
    }
}




void wrapper(double *volumeVect, double *ang0Vect,
                double *basFricAngVect, double *intAngVect, 
                double *volBasVect, double *y,
                double *circterms, double *parest, int npts,
                double *combo){
                
               
    printf("Wrapping up parameters ...\n");            
	int i;
	combo[0] = (double) npts;
	
	int n; // just to keep the indices below cleaner
	
	for(i = 1; i<= npts; i++){
		n = i-1;
		
		combo[i] = volumeVect[n];
		combo[i+(1*npts)] = ang0Vect[n];
		combo[i+(2*npts)] = basFricAngVect[n];
		combo[i+(3*npts)] = intAngVect[n];
		combo[i+(4*npts)] = volBasVect[n];
		combo[i+(5*npts)] = y[n];
	}
	//printf("Here");
	int indtemp = (npts * 6) + 1;
	n = 0;
	for (i = indtemp; i< (indtemp+5); i++){
		combo[i] = circterms[n];
		n++;
	}
	indtemp = i;
	n = 0;
	for(i = indtemp; i<(indtemp+4); i++){
		combo[i] = parest[n]; 
		n++;
	}
	
	printf("Done wrapping up parameters\n");
}








double prTimesLhood(double *volumeVect, double *ang0Vect,
                double *basFricAngVect, double *intAngVect, 
                double *volBasVect, double *y, double *thetas,
                double *circterms, double *parest, int npts) 
{  

	vect_Exp(thetas, thetas, 3);
	//printf("begining prTimesLhood: ");
	//vect_Print(thetas, 3, 4);
    double theta1 = thetas[0], theta2 = thetas[1];
    double theta3 = thetas[2];

    int ndp = npts;
    double alpha = 1.9, alphat = 2;
    double beta_v = parest[0], m_v = parest[1], m_bf = parest[2];
    double m_vbf = parest[3];

    double **C = create2Darray(ndp, ndp, 0.0);
    double **dCdbeta1 = create2Darray(ndp, ndp, 0.0); 
    double **dCdbeta2 = create2Darray(ndp, ndp, 0.0);
    double **dCdbeta3 = create2Darray(ndp, ndp, 0.0);

    double cn[5]; // !! could change depending on #cols of design points !!
    int i, j;
    for(i = 0; i<5; i++){
        cn[i] = exp((-1* pow(circterms[i],2))/(4*theta2))/sqrt(PI * theta2);
    }

    double c0 = 1/(2*sqrt(PI * theta2));
    //printf("%3.4f\n", c0);

    double dcn[5]; // !! could change depending on #cols of design points !!

    for (i = 0; i<5; i++){
        dcn[i] = (-1/(2*theta2 * sqrt(PI*theta2))) *
                (1- (pow(circterms[i], 2)/(2*theta2))) *
                exp(-1*pow(circterms[i], 2)/(4*theta2));
    }

    //vect_Print(dcn, 5);
    double dc0 = -1/(2*theta2 * sqrt(PI * theta2)) * 0.5;
    //printf("%3.4f\n", dc0);

    int jj, kk;
    double temp1[5], temp2[5];
    double corrangle;
    double dcorrangle;
	double temp1Sum;
	double temp2Sum;
	
    for (jj = 0; jj<ndp; jj++){
        for (kk= 0; kk<ndp; kk++){
            /* Calculate the correlation due to angle*/
            //printf("cn: "); vect_Print(cn, 5);
            //printf("circterms: "); vect_Print(circterms, 5);
            vect_MultN(temp1, temp1, 0, 5);
            vect_MultN(temp1, temp1, 0, 5);
            //printf("jj & kk: %3.4f  %3.4f\n", ang0Vect[jj], ang0Vect[kk]);
            for (i = 0; i<5; i++){
                temp1[i] = cn[i] * cos(circterms[i] *
                            (ang0Vect[jj] - ang0Vect[kk]));
                temp2[i] = dcn[i] * cos(circterms[i] *
                            (ang0Vect[jj] - ang0Vect[kk]));
            }

        //printf("Came back ===> ");vect_Print(temp1, 5);
        temp1Sum = vect_Sum(temp1, 5);
        temp2Sum = vect_Sum(temp2, 5);
        //printf("temp1Sum = %3.4f\n", temp1Sum);
        corrangle = c0 + temp1Sum;
        dcorrangle = dc0 + temp2Sum;

        C[jj][kk] = exp(-1*theta1 *pow(abs_double(volumeVect[jj] -
                     volumeVect[kk]), alpha)) * exp(-1*theta3 *
                     pow(abs_double(basFricAngVect[jj] -
                     basFricAngVect[kk]), alpha)) * corrangle;

		dCdbeta1[jj][kk] = pow_double(-1 * abs_double(volumeVect[jj] -
                           volumeVect[kk]), alpha) * C[jj][kk];

		dCdbeta2[jj][kk] = C[jj][kk] * dcorrangle /corrangle;
		dCdbeta3[jj][kk] = pow_double(-1 * abs_double(basFricAngVect[jj] -
                           basFricAngVect[kk]),alpha) * C[jj][kk];
		
		
        }
    }

 


    double *XX = (double*) malloc(sizeof(double) * npts);
    double *YY = (double*) malloc(sizeof(double) * npts);
    for (i= 0; i< npts; i++){
        XX[i] = beta_v + (m_v * volumeVect[i]) +
                (m_bf * basFricAngVect[i])+ (m_vbf*volBasVect[i]);
        YY[i] = y[i];
    }
 


    double *Rinv = (double*) malloc(sizeof(double) * (ndp*ndp));
    double *Rinv2 = (double*) malloc(sizeof(double) * (ndp*ndp));
    int Rinv_index = 0; 
    for (i = 0; i<ndp; i++){
        for(j = 0; j<ndp; j++){
            Rinv[Rinv_index] = C[j][i];
            Rinv2[Rinv_index] = C[j][i];
            Rinv_index++;
        }
    }    

    matrix_Inverse(Rinv, ndp);
    matrix_Inverse(Rinv2, ndp);    
    double *signhalf = (double*) malloc(sizeof(double) * (ndp*ndp));

    matrix_sqrt(Rinv, ndp, signhalf);
    //matrix_Print3(signhalf, ndp, ndp, 4);
    double *temponet = (double*) malloc(sizeof (double) * ndp);
    for (i = 0; i < ndp; i++){
        temponet[i] = 1;
    }
    double *onet = (double*) malloc(sizeof(double) * ndp);
    matrix_Mult('N', 'N', ndp, 1, ndp, 1.0, signhalf, ndp, temponet,
                ndp, 0.0, onet, ndp);
    free(temponet);
    

	/**     Calculating theta_hat    **/   
    /*ttemp = XX' * Rinv*/
    double *ttemp = (double*) malloc(sizeof(double) * ndp);
    matrix_Mult('T', 'N', 1, ndp, ndp, 1.0, XX, ndp, Rinv2, ndp, 
                0.0, ttemp, 1);

    /* templeft = XX' * Rinv * XX  : 1 number since XX is a column */
	double templeft;
    matrix_Mult('N', 'N', 1, 1, ndp, 1.0, ttemp, 1, XX, ndp, 0.0, 
                &templeft, 1);
    /* tempmid = Rinv * YY*/
    double *tempmid = (double*) malloc(sizeof(double) * ndp);
    matrix_Mult('N', 'N', ndp, 1, ndp, 1.0, Rinv2, ndp, XX, ndp, 
                0.0, tempmid, ndp); 

    /* tempright= (Rinv * YY)' * YY*/
    double tempright;
    matrix_Mult('T', 'N', 1, 1, ndp, 1.0, tempmid, ndp, YY, ndp, 0.0,
    			&tempright, 1);
     
	/* theta_hat = inv(XX' * Rinv * XX) * (Rinv * XX)' * YY */
	double theta_hat = (tempright)/(templeft);
 
    
    //printf("Theta_hat= %3.4f\n", theta_hat);
	
	// ---------------------------------------------------- //	



	/**             Calculating bweights                   **/
    double *bweights = (double*) malloc(sizeof(double) * ndp);
    double *bwtemp = (double*) malloc(sizeof(double) * ndp);

    for (i=0; i<ndp; i++){
        bwtemp[i] = YY[i] - (XX[i] * theta_hat);
    }

    matrix_Mult('N', 'N', ndp, 1, ndp, 1.0, Rinv2, ndp,bwtemp, 
                ndp, 0.0, bweights, ndp); 

	// ----------------------------------------------------- //
	
	// bw2 = Rinv * XX = C\XX
    double *bw2 = (double*) malloc(sizeof(double) * ndp);
    matrix_Mult('N', 'N', ndp, 1, ndp, 1.0, Rinv2, ndp,XX, 
                ndp, 0.0, bw2, ndp); 

    double ssqrd;
    matrix_Mult('T', 'N', 1, 1, ndp, 1.0, bwtemp, ndp, bweights,
                 ndp, 0.0, &ssqrd, 1);
    
    //printf("%3.4f\n", ssqrd);


	
	/**                 Calculating L                          **/    
    double *Cvect = (double*) malloc(sizeof(double) * (ndp*ndp));   
    int CvectIndex = 0;
    for (i = 0; i<ndp; i++){
        for(j = 0; j<ndp; j++){
            Cvect[CvectIndex] = C[j][i];
            CvectIndex++;
        }
    }
    
    double detC; 
    matrix_det(Cvect,ndp, ndp, &detC);
    //printf("determinant = %E\n", detC );
	double L = (pow_double(detC, -0.5))* (pow_double(abs(templeft), -0.5)) * 
   				(pow_double(ssqrd, (-(ndp-3)/2)));
    //printf("%E \n", L);
	// --------------------------------------------------------- //


   /**                  Calculating tempmat                    **/
	
	// XX_XXp = XX' * XX
    double *XX_XXp = (double*) malloc(sizeof(double)*(ndp*ndp));
    matrix_Mult('N','T', ndp, ndp, 1, 1.0, XX, ndp, XX, ndp, 0.0,
                XX_XXp, ndp);
    
    // ttempmat = -XX * inv(XX' * Rinv * XX) * XX'
    double *ttempmat = (double*) malloc(sizeof(double)*(ndp*ndp));
    for(i=0; i<(ndp*ndp); i++){
        ttempmat[i] = 1;
    }
    matrix_Mult('N','N', ndp, ndp, ndp, (-1.0/(templeft)), XX_XXp,
                ndp, Rinv2, ndp, 1.0, ttempmat, ndp);
     
    //matrix_PrintEnd(ttempmat, ndp*ndp, 10);

    /* tempmat = Rinv * ((ones(ndp)) - XX*inv(XX'*Rinv*XX)*XX'*Rinv) */
    double *tempmat = (double*) malloc(sizeof(double)*(ndp*ndp));
    matrix_Mult('N','N', ndp, ndp, ndp, 1.0, Rinv2, ndp, ttempmat, 
                ndp, 0.0, tempmat, ndp);
    
    //matrix_PrintEnd(tempmat, ndp*ndp, 10);
    
    // --------------------------------------------------------- //
    double *w1 = (double*) malloc(sizeof(double)*(ndp*ndp));
    double *w2 = (double*) malloc(sizeof(double)*(ndp*ndp));
    double *w3 = (double*) malloc(sizeof(double)*(ndp*ndp));
    double *dCdbeta1Vect = (double*) malloc(sizeof(double)*(ndp*ndp));
    double *dCdbeta2Vect = (double*) malloc(sizeof(double)*(ndp*ndp));
    double *dCdbeta3Vect = (double*) malloc(sizeof(double)*(ndp*ndp));
    int vectInd = 0;
    for (i = 0; i<ndp; i++){
        for(j = 0; j<ndp; j++){
            dCdbeta1Vect[vectInd] = dCdbeta1[j][i];
            dCdbeta2Vect[vectInd] = dCdbeta2[j][i];
            dCdbeta3Vect[vectInd] = dCdbeta3[j][i];
            vectInd++;
        }
    }
    
    matrix_Mult('N','N', ndp, ndp, ndp, 1.0, dCdbeta1Vect, ndp, tempmat,
                ndp, 0.0, w1, ndp);
    matrix_Mult('N','N', ndp, ndp, ndp, 1.0, dCdbeta2Vect, ndp, tempmat,
                ndp, 0.0, w2, ndp);
    matrix_Mult('N','N', ndp, ndp, ndp, 1.0, dCdbeta3Vect, ndp, tempmat,
                ndp, 0.0, w3, ndp);

    //matrix_PrintEnd(w3, ndp*ndp, 5);
    double *Istar = (double*) malloc(sizeof(double)*(4*4));
    Istar[0] = ndp-3;
    matrix_Trace(w1, &Istar[4], ndp);
    Istar[1] = Istar[4];
    matrix_Trace(w2, &Istar[8], ndp);
    Istar[2] = Istar[8];
    matrix_Trace(w3, &Istar[12], ndp);
    Istar[3] = Istar[12];
    
    double *w1w1 = (double*) malloc(sizeof(double) * (ndp*ndp));
    matrix_Mult('N','N', ndp, ndp, ndp, 1.0, w1,ndp, w1, ndp, 0.0,
                w1w1, ndp);
    matrix_Trace(w1w1, &Istar[5], ndp);
    

    double *w2w2 = (double*) malloc(sizeof(double) * (ndp*ndp));
    matrix_Mult('N','N', ndp, ndp, ndp, 1.0, w2,ndp, w2, ndp, 0.0,
                w2w2, ndp);
    matrix_Trace(w2w2, &Istar[10], ndp);
    
     
    double *w3w3 = (double*) malloc(sizeof(double) * (ndp*ndp));
    matrix_Mult('N','N', ndp, ndp, ndp, 1.0, w3,ndp, w3, ndp, 0.0,
                w3w3, ndp);
    matrix_Trace(w3w3, &Istar[15], ndp);
    
    double *w1w2 = (double*) malloc(sizeof(double) * (ndp*ndp));
    matrix_Mult('N','N', ndp, ndp, ndp, 1.0, w1,ndp, w2, ndp, 0.0,
                w1w2, ndp);
    matrix_Trace(w1w2, &Istar[9], ndp);

    double *w1w3 = (double*) malloc(sizeof(double) * (ndp*ndp));
    matrix_Mult('N','N', ndp, ndp, ndp, 1.0, w1,ndp, w3, ndp, 0.0,
                w1w3, ndp);
    matrix_Trace(w1w3, &Istar[13], ndp);
    
    Istar[6] = Istar[9];
    
    Istar[7] = Istar[13];
    
    double *w2w3 = (double*) malloc(sizeof(double) * (ndp*ndp));
    matrix_Mult('N','N', ndp, ndp, ndp, 1.0, w2,ndp, w3, ndp, 0.0,
                w2w3, ndp);
    matrix_Trace(w2w3, &Istar[14], ndp);
    Istar[11] = Istar[14];  

    //matrix_Print3(Istar, 4, 4, 4);
    
    double det_Istar;
    matrix_det(Istar, 4,4, &det_Istar);
    double thetas_prod;
    vect_Prod(thetas, &thetas_prod, 3);
    double output;    
    output = (-1)* L* sqrt(abs(det_Istar)) * thetas_prod;

   // printf("This output = %E\n", (output));
    
    if (isnan(output) || isinf(output)){
        output = 0.0;
    }

    return output;

}


double wrprTimesLhood(double *thetas, double *wrparams){

	int npts = (int) wrparams[0];
	int i;
	
	double *vv = (double*) malloc(sizeof(double) * npts);
	double *av = (double*) malloc(sizeof(double) * npts);
	double *bfav = (double*) malloc(sizeof(double) * npts);
	double *iav = (double*) malloc(sizeof(double) * npts);
	double *vbv = (double*) malloc(sizeof(double) * npts);
	double *y = (double*) malloc(sizeof(double) * npts);
	
	double *circ = (double*) malloc(sizeof(double) * 5);
	double *par = (double*) malloc(sizeof(double) * 5);
	
	int n;
	
	for(i = 1; i<= npts; i++){
		n = i-1;
		vv[n] = wrparams[i];
		av[n] = wrparams[i+(1*npts)];
		bfav[n] = wrparams[i+(2*npts)];
		iav[n] = wrparams[i+(3*npts)];
		vbv[n] = wrparams[i+(4*npts)];
		y[n] = wrparams[i+(5*npts)];
	}
	
	int indtemp = (npts * 6) + 1;
	n = 0;
	for (i = indtemp; i< (indtemp+5); i++){
		circ[n] = wrparams[i];
		n++;
	}
	
	indtemp = i;
	n = 0;
	for(i = indtemp; i<(indtemp+4); i++){
		par[n] = wrparams[i]; 
		n++;
	}
	
	//vect_Print(par, 4, 4);
	
	double result;
	
	result = prTimesLhood(vv, av, bfav, iav, vbv, y, thetas, circ, 
							par, npts);
	
	free(vv); free(av); free(bfav); free(iav); free(vbv);
	free(y); free(circ); free(par);
	
	
	return result;
	
	
}

double gsl_wrprTimesLhood(const gsl_vector *thetasvect, void * wrparams){
/* Evaluate the function wrprTimesLhood using a gsl_vector instead of a 
double* 
*/
	double *thetas = (double*) malloc(sizeof(double) * 3);
	int i;
	for(i = 0; i<3; i++){
		thetas[i] = gsl_vector_get(thetasvect, i);
	}
	
	double output = wrprTimesLhood(thetas, wrparams);
	free(thetas);
	return output;
	
}





void posteriorMode(double *thetas, double *wrparams){
	
	//int i; 
	//for(i =0;i<3; i++){
	//	printf("theta[%d]= %3.4f\n", i, thetas[i]);
	//}
	//vect_Print(thetas, 3, 4);
	//vect_Print(wrparams, (207*6)+10, 4);
	//clock_t tic = clock();
	
	const gsl_multimin_fminimizer_type *T =
			gsl_multimin_fminimizer_nmsimplex;

	gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

	int iter = 0; int max_iter = 100;
    int status;
    double size;
	int n = 3; //length of variable output. in our case: 3.
	/*STARTING POINT*/
	x = gsl_vector_alloc(n);
	gsl_vector_set(x, 0, thetas[0]);
	gsl_vector_set(x, 1, thetas[1]);
	gsl_vector_set(x, 2, thetas[2]);
	
	/*SET INITIAL STEP SIZE TO 1*/
	ss = gsl_vector_alloc(n);
	gsl_vector_set_all(ss, 1.0);
	
	/*INITIALIZE METHOD AND ITERATE*/
	minex_func.n = n;
	minex_func.f = gsl_wrprTimesLhood;
	minex_func.params = wrparams;
	
	double tot = gsl_wrprTimesLhood(x, wrparams);
	//printf("tot = %E\n", tot);
	
	s = gsl_multimin_fminimizer_alloc(T,n);
	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
	
	
	double tol = 1e-3; // tolerance.
	
	do{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		
		if (status){
			printf("break\n");
			break;
		}
		
		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, tol);
		
		if (status == GSL_SUCCESS){
    //        printf("Converged to minimum at\n");
        }
		
		// Print every 5 iterations
		if ((iter % 5) == 0){
        	printf ("%5d %3.4f %3.4f %3.4f f() = %E size = %E\n",iter,
                gsl_vector_get (s->x, 0),gsl_vector_get (s->x, 1),
                gsl_vector_get (s->x, 2), s->fval, size);
			}
        //printf("check = %d  && %d\n", (status == GSL_CONTINUE), iter<100);
    }
    while (status == GSL_CONTINUE && iter < max_iter);



	int i;
	for (i = 0; i<3; i++){
		thetas[i] = gsl_vector_get(s->x, i);
	}
	//clock_t toc = clock();
	//printf("Minimization process took %f secs\n",
		//	(double)(toc - tic) / CLOCKS_PER_SEC);
	//printf("\nLeaving posteriorMode\n");

}
















