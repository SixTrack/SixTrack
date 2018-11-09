#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//----------------------------------------------------------
//read matrix : cin>> in c++ , scanf() in c
int cin(float a[100][100]){
	int i,j,n;
	printf("\n Enter Length Of Matrix N*N : ");
	scanf("%d",&n);
	printf("\n--------------------------\n");
	for(i=0;i<n;i++)
		for(j=0;j<n;j++){
			printf(" Matrix[%d][%d] : ",i+1,j+1);
			scanf("%f",&a[i][j]);
		}
	printf("\n----------------------------------------------------\n");
return n;
}

//-----------------------------------------------------
// show matrix : cout<< in c++ , printf() in c
void cout(float a[100][100],int n,int show){
	int i,j;
	if(show == 1)
		for(i=0;i < n;i++){
			for(j=0;j < n;j++)
				printf(" %.2f \t",a[i][j]);
			printf("\n");
		}
	else if(show == 2){
		printf("\n\n The Inverse Of Matrix Is : \n\n");
		for (i=0;i<n;i++){
			for (j=0;j<n;j++)
				printf(" %.4f \t",a[i][j]);
			printf("\n");
		}
	}
}

//---------------------------------------------------
//	calculate minor of matrix OR build new matrix : k-had = minor
void minora(float b[100][100],float a[100][100],int i,int n){
	int j,l,h=0,k=0;
	for(l=1;l<n;l++)
		for( j=0;j<n;j++){
			if(j == i)
				continue;
			b[h][k] = a[l][j];
			k++;
			if(k == (n-1)){
				h++;
				k=0;
			}
		}
}// end function

//---------------------------------------------------
//	calculate determinte of matrix
float det(float a[100][100],int n){
	int i;
	float b[100][100],sum=0;
	if (n == 1)
return a[0][0];
	else if(n == 2)
return (a[0][0]*a[1][1]-a[0][1]*a[1][0]);
	else
		for(i=0;i<n;i++){
			minora(b,a,i,n);	// read function
			sum = (float) (sum+a[0][i]*pow(-1,i)*det(b,(n-1)));	// read function	// sum = determinte matrix
		}
return sum;
}// end function

//---------------------------------------------------
//	calculate transpose of matrix
void transpose(float c[100][100],float d[100][100],float n,float det){
	int i,j;
	float b[100][100];
	for (i=0;i<n;i++)
		for (j=0;j<n;j++)
			b[i][j] = c[j][i];
	for (i=0;i<n;i++)
		for (j=0;j<n;j++)
			d[i][j] = b[i][j]/det;	// array d[][] = inverse matrix
}// end function

//---------------------------------------------------
//	calculate cofactor of matrix
void cofactor(float a[100][100],float d[100][100],float n,float determinte){
	float b[100][100],c[100][100];
	int l,h,m,k,i,j;
	for (h=0;h<n;h++)
		for (l=0;l<n;l++){
			m=0;
			k=0;
			for (i=0;i<n;i++)
				for (j=0;j<n;j++)
					if (i != h && j != l){
						b[m][k]=a[i][j];
						if (k<(n-2))
							k++;
						else{
							k=0;
							m++;
						}
					}
			c[h][l] = pow(-1,(h+l))*det(b,(n-1));	// c = cofactor Matrix
		}
	transpose(c,d,n,determinte);	// read function
}// end function

//---------------------------------------------------
//	calculate inverse of matrix
void inverse(float a[100][100],float d[100][100],int n,float det){
	if(det == 0)
		printf("\nInverse of Entered Matrix is not possible\n");
	else if(n == 1)
		d[0][0] = 1;
	else
		cofactor(a,d,n,det);	// read function
}// end function

//---------------------------------------------------
//main fuction exe
int main(void){
	int i,j,n;
	float a[100][100],d[100][100],deter;
	printf("\n C Program To Find Inverse Of Matrix\n\n"); 
	n = cin(a);	// read function
	int print_matrix = 1;
	cout(a,n,print_matrix);	// read function
	deter = (float) det(a,n);	// read function
		printf("----------------------------------------------------\n");
		printf("\n\n Determinant of the Matrix : %.4f ",deter);
		printf("\n\n-----------------------\n");
	inverse(a,d,n,deter);	// read function
	int print_inverse = 2;
	cout(d,n,print_inverse);	// read function
		printf("\n\n==============================* THE END *==============================\n");
		printf("\n		**** Thanks For Using The Program!!! ****\n");
	//getc();
    return 0;
}// end main
