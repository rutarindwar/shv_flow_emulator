#ifndef MATRIXOPERATIONS_H_INCLUDED
#define MATRIXOPERATIONS_H_INCLUDED

#define PI 3.14159265358979323846

void vect_Exp(double *lst, double *result, int len);
void vect_Print(double *lst, int len, int precision);
void vect_PrintINT(int *lst, int len);
void vect_PrintPart(double *lst, int a, int b, int precision);

void vect_MultN(double *lst, double *result, double k, int len);
void matrix_Print1(double* mat, int nrows, int ncols);
double ***create3Darray(int nrows, int ncols, int nplanes, double val);
double ***create3Darray2(int nplanes, int nrows, int ncols, double val);


double **create2Darray(int nrows, int ncols, double val);
void free2Darray(double **mat, int nrows, int ncols);
void matrix_Print2 (double* mat,int nrows, int ncols, int precision);
void matrix_MultN(double* mat, double k, int nrows, int ncols);


void vect_Power(double *lst, double *result, int k, int len);
double vect_Sum(double *lst, int len);

void matrix_PrintPart(double* mat, int nrows, int ncols, int ri,
                      int rj, int ci, int cj, int precision);
void matrix_PrintPart2(double **mat, int c0, int c1, int r0, int r1,
						int precision);
void matrix_PrintPart3(double **mat, int c0, int c1, int r0, int r1);
void matrix_To_Vector(double* mat, double lst[], int nrows, 
                      int ncols, int len);
                      
                      

void matrix_Print3(double* mat, int nrows, int ncols, int precision);
void matrix_Print4(double **mat, int nrows, int ncols, int precision);
void matrix_Inverse(double* mat, int size);
void matrix_sqrt(double* mat, int size, double* sqmat);

void matrix_Mult(char tra, char trb, int m, int n, int k,
                  double al, double *mata, int lda,
                  double *matb, int ldb, double bet,
                  double *matc, int ldc);


void matrix_det(double *mat, int nrows, int ncols, double *det);

void matrix_PrintEnd(double *mat, double size, int len);
void matrix_Trace(double *mat, double *trace, int size);
void vect_Prod(double *lst, double *prod, int len);

void vect_LinFit(double *x, double *y, int len, double *coef);
void matrix_solveAxb(double *a, double *b, int size);
void vect_Linspace(double *lst, double a, double b, int len);
double vect_Max(double *lst, int len);
double vect_Min(double *lst, int len);
int compare (const void *a, const void *b);
int compareINT(const void *a, const void *b);
void vect_Sort(double *lst, int *inds, int len);
void vect_Sort2(double *lst, double *sortedLst, int *inds, int len);
void matrix_Copy3Darray(double ***source, double ***dest, int nr, int nc, 
						int np);


void matrix_Transpose(double **mat, int nrows, int ncols);
void matrix_ZeroOut(double **mat, int nrows, int ncols);
int lookUpIndex(double x, double *foo, int len);
void vect2DInterp(double *xvect, double *yvect, double **zgrid, double *xi,
					double *yi, double *zi, int len);
					

double matrix_Min(double **mat, int nrows, int ncols);

void matrix_PrintToFileDOUBLE(double **mat, FILE *fout, int nrows, 
								int ncols, int precision);
void vect_PrintToFileDOUBLE(double *vect, FILE *fout, int nrows, 
								int precision);

void vect_checkMonot(double *vect, int len);

int findIndex(double x, double *xlist, int len, int *index);

void vect_Interpolation(double *xvect, double *yvect, double *xi, double *yi,
							int xvect_len, int xi_len);



double matrix_Max(double **mat, int nrows, int ncols);


double matrix_Sum(double **mat, int nrows, int ncols);


void vect2DInterp2(double *xvect, double *yvect, double **zgrid, double *xi,
					double *yi, double **zi, int xvlen, int yvlen, 
					int xilen, int yilen);


void matrix_readFromFile(FILE *fp, double **mat, int nrows, int ncols);
 









//int compare (const void *const first, const void *const second);
#endif // MATRIXOPERATIONS_H_INCLUDED
