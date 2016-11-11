#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>


#include "matrixOperations.h"
#include "readInZgrid.h"
#include "f2c.h"
#define MAXNAME 1000
#define MAXLINEBUFFER 5000


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

void readVectFromFileDOUBLE(FILE *fp, double *vect, int size){
    int i;
    for(i = 0; i <size; i++){
        fscanf(fp, "%lf", &vect[i]);
    }

}






int main(int argc, char *argv[]){

    char fnameIndListN[MAXNAME];
    char fnameIndListW[MAXNAME];
    strcpy(fnameIndListN, "./newrunsdir_BV_100m/indlistnorth.txt");
    strcpy(fnameIndListW, "./newrunsdir_BV_100m/indlistwest.txt");
    FILE *findlstN;
    FILE *findlstW;

    findlstN = fopen(fnameIndListN, "r");
    findlstW = fopen(fnameIndListW, "r");
    int nbrIndLst = linesInFile(findlstN);
    //printf("nbrIndlist = %d\n", nbrIndLst);
    int *indN = (int*) malloc(sizeof(int) * nbrIndLst);
    int *indW = (int*) malloc(sizeof(int) * nbrIndLst);
    readVectFromFileINT(findlstN, indN, nbrIndLst);
    readVectFromFileINT(findlstW, indW, nbrIndLst);
    fclose(findlstN);
    fclose(findlstW);

   // vect_PrintINT(indN, nbrIndLst);
    int kmax = nbrIndLst;
    
    char dirname[MAXNAME];
    strcpy(dirname, "./newrunsdir");
    char filename[MAXNAME];

    int i; int j;
    int k;

    /*File-dependent variables*/
    FILE *fheight;
    snprintf(filename, MAXNAME, "%s/heights_%06d.txt", dirname,1);
    fheight = fopen(filename, "r");
    int hvlen = linesInFile(fheight);
    fclose(fheight);

    double *y = (double*) malloc(sizeof(double) * hvlen);
    double *H = (double*) malloc(sizeof(double) * hvlen);
    
    /* Reading in danilo_qs.txt*/
    int nc; int nr;
    FILE *fdan = fopen("danilo_qs.txt", "r");
    nc = colsInFile(fdan);
    nr = linesInFile(fdan);
    //printf("danilo_qs has %d rows and %d cols\n", nr, nc);
    double **danqs = create2Darray(nr, nc, 0.0);
    matrix_readFromFile(fdan,danqs, nr, nc); 
    fclose(fdan);
    //matrix_PrintPart2(danqs, nr-3, nr-1, 0, 8, 4);

    /* Reading in uncertain_input_list*/
    snprintf(filename, MAXNAME, "%s/uncertain_input_list.txt", dirname);
    //puts(filename);
    FILE *fuil = fopen(filename, "r");
    nc = colsInFile(fuil);
    nr = linesInFile(fuil);
    //printf("uncertain_input_list has %d rows and %d cols\n", nr, nc);
    double **uilist = create2Darray(nr, nc, 0.0);
    matrix_readFromFile(fuil, uilist, nr, nc);
    fclose(fuil);
    //matrix_PrintPart2(uilist, nr-3, nr-1, 0, 3, 8);
   
    
    double *uilist1 = (double*) malloc(sizeof(double) * nr);
    for(i = 0; i<nr; i++){
        uilist1[i] = uilist[i][0];
    }
    //vect_Print(uilist1, nr, 4);
    
    double *tvols = (double*) malloc(sizeof(double) * nr);
    int *inds = (int*) malloc(sizeof(int) * nr);
    vect_Sort2(uilist1,tvols, inds, nr);
    //vect_PrintPart(tvols, 0, 2, 4);
    //vect_PrintINT(inds, nr);
    
    int indslen = nr;
    int *indsavetemp = (int*) malloc(sizeof(int) * indslen);
    
    int indsavelen = 0;

    //printf("indslen = %d \n", indslen);
    for(j = 0; j<indslen; j++){
        if ((uilist[inds[j]][2] > danqs[j][5]) &&
            (uilist[inds[j]][2] < danqs[j][6])){

            indsavetemp[indsavelen] = j;
            indsavelen++;
        }
    }

    printf("indsavelen = %d\n", indsavelen);
    int *indsave = (int*) malloc(sizeof(int) * indsavelen);

    for(i = 0; i<indsavelen; i++){
        indsave[i] = indsavetemp[i];
    }
    //int *indsave = (int*) realloc(indsavetemp, indsavelen*sizeof(int));
    //vect_PrintINT(indsave, indsavelen);
    free(indsavetemp);
    
    /*----------------------------------------------------------------
     * Arrays angs, bangs, iangs, nbangs, ... aren't used in the code.
     * Implement them later as needed.
     * --------------------------------------------------------------*/
    
    double *tangs = (double*) malloc(sizeof(double) * indslen);
    double *tbangs = (double*) malloc(sizeof(double) * indslen);
    double *tiangs = (double*) malloc(sizeof(double) * indslen);

    for(i = 0; i< indslen; i++){
        tangs[i] = uilist[inds[i]][1];
        tbangs[i] = uilist[inds[i]][2];
        tiangs[i] = uilist[inds[i]][3];
    }
    free2Darray(uilist, indslen, nc);
    //vect_PrintPart(tiangs, 0, 9, 4);
    
    double *angs = (double*) malloc(sizeof(double) * indslen);
    double *bangs = (double*) malloc(sizeof(double) * indslen);
    double *iangs = (double*) malloc(sizeof(double) * indslen);
    double *vols = (double*) malloc(sizeof(double) * indslen);

    for(i = 0; i<indsavelen; i++){
        angs[i] = tangs[indsave[i]];
        bangs[i] = tbangs[indsave[i]];
        iangs[i] = tiangs[indsave[i]];
        vols[i] = tvols[indsave[i]];
    }

   /* Normalizing input variables*/ 
    double *nangs = (double*) malloc(sizeof(double) * indsavelen);
    double *nvols = (double*) malloc(sizeof(double) * indsavelen);
    double *nbangs = (double*) malloc(sizeof(double) * indsavelen);

    for(i = 0; i<indsavelen; i++){
        nangs[i] = angs[i] / 360.0;
        nvols[i] = (vols[i] - vect_Min(vols, indsavelen))/
                    (vect_Max(vols, indsavelen)-
                                vect_Min(vols, indsavelen));
    
        nbangs[i] = log(tan(bangs[i] * PI/180));
    }

    double maxtmp = vect_Max(nbangs, indsavelen);
    double mintmp = vect_Min(nbangs, indsavelen);

    for(i = 0; i<indsavelen; i++){
        nbangs[i] = abs(nbangs[i]/(maxtmp - mintmp));
    }
    
    mintmp = vect_Min(nbangs, indsavelen);
    for(i = 0; i<indsavelen; i++){
        nbangs[i] = nbangs[i] - mintmp;
    }
    //vect_PrintPart(nbangs, 0, 9, 4);


    double *th = (double*) malloc(sizeof(double) * hvlen);
    double *hindsave = (double*) malloc(sizeof(double) * indsavelen);
    
    int posinds;
    double abound = 0.10;
    double bbound = 2.0;
    int designYrows = 0;
    int designXrows = 0;
    int designXcols = 5;
    
    double **designX = create2Darray(indsavelen, designXcols, 0.0);
    double *designY = (double*) malloc(sizeof(double) * indsavelen);
    double designYtemp;

    int *smallinds = (int*) malloc(sizeof(int) * indsavelen);
    int smallindslen;
    int *zeroinds = (int*) malloc(sizeof(int) * indsavelen);
    int zeroindslen;
    

    double dist2z; double dist2zp1; double dist2zm1;

    int minindex;
    double tmpmin;

    char outdirname[MAXNAME];
    strcpy(outdirname, "./subdesignsIO");

    FILE *fpdesignY;
    FILE *fpdesignX;

    k = atoi(argv[1]) + 1;
   // for (k = 1; k <= 1953; k++){
        
    snprintf(filename, MAXNAME, "%s/heights_%06d.txt", dirname,k);
    puts(filename);
    fheight = fopen(filename, "r");
    readVectFromFileDOUBLE(fheight, y, hvlen);
    //vect_Print(y, hvlen, 4); 
    fclose(fheight);
    for(i =0; i<hvlen; i++){
        H[i] = log(1 + y[i]);
    }
   // vect_PrintPart(H, 65, 67, 4);
    for(i = 0; i< hvlen; i++){
        th[i] = H[inds[i]];
    }

    //vect_PrintPart(th, 1567, 1570, 4);
   
    for(i = 0; i<indsavelen; i++){
        hindsave[i] = th[indsave[i]];
    }
    
    //vect_PrintPart(hindsave, indsavelen-3, indsavelen-1, 4);
        
    posinds = 0;
    for(i = 0; i<indsavelen; i++){
        if (hindsave[i] > log(1.10)){
            posinds++;
        }
    }
    printf("posinds = %d\n", posinds);

    if((posinds == 0) || (posinds < 15)){
        designYrows = 1;
        designXrows = 1;
    }
    //vect_Print(designY, designYrows, 4);
    //matrix_Print4(designX,designXrows, designXcols, 4);    
     
    else{
        //vect_Print(hindsave, indsavelen, 4);
        smallindslen = 0;
        zeroindslen = 0;
        for(i = 0; i<indsavelen; i++){
            if((hindsave[i] > log(1+abound)) &&
                (hindsave[i] < log(bbound))){
                smallinds[smallindslen] = i;
                smallindslen++;
            }
            if (hindsave[i] < log(1 + abound)){
                zeroinds[zeroindslen] = i;
                zeroindslen++;
            }
        }
        //vect_PrintINT(zeroinds, zeroindslen);
        double *smallestdist2z = (double*) malloc(sizeof(double) * 
                                                    zeroindslen);
                                                                
        int *zerokeep = (int*) malloc(sizeof(int*)* 
                                           smallindslen);
        int zerokeeplen = 0;

        if (zeroindslen >= 1){
           zerokeeplen = smallindslen;
            for(j = 0; j<smallindslen; j++){
                tmpmin = 100.0; // VERY ARBITRARY VALUE.
                for(i = 0; i<zeroindslen; i++){
                
                    dist2z= pow(nbangs[zeroinds[i]] - 
                                nbangs[smallinds[j]], 2) +
                            pow(nangs[zeroinds[i]] - 
                                nangs[smallinds[j]], 2) +
                            pow(nvols[zeroinds[i]] - 
                                nvols[smallinds[j]], 2);

                    dist2zp1= pow(nbangs[zeroinds[i]] - 
                                nbangs[smallinds[j]], 2) +
                            pow(nangs[zeroinds[i]] + 1 - 
                                nangs[smallinds[j]], 2) +
                            pow(nvols[zeroinds[i]] - 
                                nvols[smallinds[j]], 2);
                    
                    dist2zm1= pow(nbangs[zeroinds[i]] - 
                                nbangs[smallinds[j]], 2) +
                            pow(nangs[zeroinds[i]] - 1 -
                                nangs[smallinds[j]], 2) +
                            pow(nvols[zeroinds[i]] - 
                                nvols[smallinds[j]], 2);
                    
                    smallestdist2z[i] = min(min(dist2z, dist2zp1),
                                            dist2zm1);
                    
                    if (smallestdist2z[i] < tmpmin){
                        minindex = i;
                        tmpmin = smallestdist2z[i];
                    }
                }

                zerokeep[j] = minindex;
            }
        }
        //vect_PrintPart(smallestdist2z, 0, 2, 4);
        //vect_PrintINT(zerokeep, smallindslen);
       // printf("zerokeeplen = %d\n", zerokeeplen);
        
        int *zinds;
        int zindslen = 0;
        if (zerokeeplen >= 1){
            qsort(zerokeep, zerokeeplen, sizeof(int), compareINT);
            zinds = (int*) malloc(sizeof(int) * zerokeeplen);
            zinds[0] = zeroinds[zerokeep[0]];
            zindslen = 1;
            //vect_PrintINT(zerokeep, zerokeeplen);
            for(j = 1; j<zerokeeplen; j++){
                if(zinds[zindslen-1] != zeroinds[zerokeep[j]]){
                    zinds[zindslen] = zeroinds[zerokeep[j]];
                    zindslen++;
                }
                qsort(zinds, zindslen, sizeof(int), compareINT);
            }
        }
        //vect_PrintINT(zinds, zindslen);
        

        int *pinds;
        int pindslen;
        if(zerokeeplen >= 1){
            pindslen = smallindslen + zindslen;
            pinds= (int*) malloc(sizeof(int)*pindslen);
            for(i = 0; i<zindslen; i++){
                pinds[i] = zinds[i];
            }
            for(i=0; i<smallindslen; i++){
                pinds[i+zindslen] = smallinds[i];
            }
        }
        if(zerokeeplen == 0){
            pindslen = smallindslen;
            pinds = (int*) malloc(sizeof(int) * pindslen);
            for(i = 0; i< pindslen; i++){
                pinds[i] = smallinds[i];
            }
        }
        //vect_PrintINT(pinds, 2);
        
        double maxv = vols[pinds[0]];
        for(i = 1; i<pindslen; i++){
            if (vols[pinds[i]] > maxv){
                maxv = vols[pinds[i]];
            }
        }
        maxv = maxv + 0.1;
        printf("maxv = %3.4f\n", maxv);
        
        int *goodindstmp = (int*) malloc(sizeof(int) * indsavelen);
        int goodindstmplen = 0;
        for(i = 0; i<indsavelen; i++){
            if((vols[i]<maxv) && (hindsave[i] > log(1 + abound))){
                goodindstmp[goodindstmplen] = i;
                goodindstmplen++;
            }
        }
        
        //printf("goodindstmplen = %d\n", goodindstmplen);
        //vect_PrintINT(goodindstmp, goodindstmplen);
        
        /* Checking if goodinds needs to be updated*/
        int goodindslen = goodindstmplen;
        if (zerokeeplen >= 1){
            goodindslen = goodindstmplen + zindslen;
        }
        
        /* Updating goodinds */
        int *goodinds = (int*) malloc(sizeof(int) * goodindslen);
        
        if (goodindslen == goodindstmplen){
            for(i = 0; i< goodindstmplen; i++){
                goodinds[i] = goodindstmp[i];
            }
        }
        if (goodindslen != goodindstmplen){
            for(i = 0; i<goodindstmplen; i++){
                goodinds[i] = goodindstmp[i];
            }
            for(i = 0; i<zindslen; i++){
                goodinds[i+goodindstmplen] = zinds[i];
            }
        }
        
        //vect_PrintINT(goodinds, goodindslen);
        printf("goodindslen = %d\n", goodindslen);

        //int np = goodindslen;
        //double fac = max(min((1/50.0) * (np - 150), 1), 0);
       // printf("fac = %3.4f\n", fac);
        

        int maxdespts = 200 ;
        if (goodindslen > maxdespts){
            printf("had more than %d design points\n", maxdespts);
            designYrows = maxdespts;
            designXrows = maxdespts;
            int tmpind;
            for(i = 0; i < maxdespts; i++){
                tmpind = rand() % (maxdespts);
                designY[i] = hindsave[goodinds[tmpind]];
                
                designX[i][0] = vols[goodinds[tmpind]];
                designX[i][1] = angs[goodinds[tmpind]];
                designX[i][2] = bangs[goodinds[tmpind]];
                designX[i][3] = iangs[goodinds[tmpind]];
                designX[i][4] = bangs[goodinds[tmpind]] * 
                                vols[goodinds[tmpind]];
            }
        }



        else{
            printf("had less than %d design points\n", maxdespts);
            designYrows = goodindslen;
            designXrows = goodindslen;
            for(i = 0; i < designYrows; i++){
                designY[i] = hindsave[goodinds[i]];
                
                designX[i][0] = vols[goodinds[i]];
                designX[i][1] = angs[goodinds[i]];
                designX[i][2] = bangs[goodinds[i]];
                designX[i][3] = iangs[goodinds[i]];
                designX[i][4] = bangs[goodinds[i]] * 
                                vols[goodinds[i]];
            }
        }

            
    } // if ..
   // printf("indsavelen= %d, goodindslen = %d \n",
   // indsavelen,designYrows);
    if (designYrows > 1){
        printf("nrows = %d\n", designYrows);
        //matrix_PrintPart2(designX, designXrows-3, designXrows-1,0,4,4);
        snprintf(filename, MAXNAME, "%s/logheightp_%06d.txt",outdirname,k);
        puts(filename);
  
        // vect_Print(designY, designYrows, 3);
        fpdesignY = fopen(filename, "w");
        vect_PrintToFileDOUBLE(designY, fpdesignY, designYrows, 4);
        fclose(fpdesignY);
        snprintf(filename, MAXNAME,"%s/designp_%06d.txt", outdirname,k);
        puts(filename);
        fpdesignX = fopen(filename, "w");
        matrix_PrintToFileDOUBLE(designX, fpdesignX, designXrows, 
                                designXcols, 4);
        fclose(fpdesignX);
    }
    else{
        printf("nrows = %d\n", designYrows);
        //matrix_PrintPart2(designX, designXrows-3, designXrows-1,0,4,4);
        snprintf(filename, MAXNAME, "%s/logheightp_%06d.txt",outdirname,k);
        puts(filename);
  
        // vect_Print(designY, designYrows, 3);
        fpdesignY = fopen(filename, "w");
        designY[0] = 0;
        vect_PrintToFileDOUBLE(designY, fpdesignY, designYrows, 4);
        fclose(fpdesignY);
        snprintf(filename, MAXNAME,"%s/designp_%06d.txt", outdirname,k);
        puts(filename);
        fpdesignX = fopen(filename, "w");
        
        designX[0][0] = 0.0; designX[0][1] = 0.0;
        designX[0][2] = 0.0; designX[0][3] = 0.0;
        designX[0][4] = 0.0;
        matrix_PrintToFileDOUBLE(designX, fpdesignX, designXrows, 
                                designXcols, 4);
        fclose(fpdesignX);

    }
    printf("_____________________________\n");



























	return 0;
}
