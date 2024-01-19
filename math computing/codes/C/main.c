#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libs/matfun.h"
#define TOLERANCE 1e-6

int  main()
{

FILE *fp; //file pointer
double **v1,**v2,**v3,**A,**B,**C,**D,**M,**F,**G,**vertices; //declaring matrices names
double AB,AC,AD,AE,BC,BD,BE,CD,CE,DE; //side lengths
int m =3, n=3; // (mxn) matrix
double l = 9; //length of a side 
double CMA,BMD,DBC,BCA;
double theta,Alpha,Beta,Beta1;

vertices = createMat(m,3);
vertices = loadMat("Vec.dat",m, 3);
printf("\n vertices\n");
printMat(vertices,3,3);

v1=createMat(m,n); //v1
v2=createMat(m,n); //v2
v3=createMat(m,n); //v3
v1=Matcol(vertices,3,0);
v2=Matcol(vertices,3,1);
v3=Matcol(vertices,3,2);

/*printf("\n v1 \n");
printMat(v1,3,1);
printf("\n v2 \n");
printMat(v2,3,1);
printf("\n v3 \n");
printMat(v3,3,1);*/

double vn1,vn2,vn3,**b1,**b11;

vn1 = norm(v1);
//printf("\n norm of v1 vn1 = %lf\n\n ",vn1);
A=Matscale(v1,3,1,1/vn1); //A= v1/vn1
printf("\n A =\n");
printMat(A,3,1);

// B

b1=tr_mul_mat(v2,A,A);// (v2.T @A )* A
b11=Matsub(v2,b1,3,1); // v2- b111
vn2=norm(b11); // norm b11
//printf("\n norm of b1 = %lf \n",vn2);

B= Matscale(b11,3,1,1/vn2);
printf("\n B\n");
printMat(B,3,1);

// C
double **v3_a,**v3_b,**c1,**c11;
v3_a=tr_mul_mat(v3,A,A);//(v3.T @ A )* A
v3_b=tr_mul_mat(v3,B,B);//(v3.T @  B)* B

c11=Matsub(v3,v3_a,3,1);
c1=Matsub(c11,v3_b,3,1);
vn3=norm(c1);
//printf("\n norm of c1 vn3 = %lf\n\n ",vn3);
C= Matscale(c1,3,1,1/vn3);

printf("\n C \n");
printMat(C,3,1);

double **ABC, **ABC_T,**I;
ABC=p_m(A,B,C,3,3); // [A B C]
printf("\n [ A B C ] \n");
printMat(ABC,3,3);

ABC_T=transposeMat(ABC,3,3);
printf("\n [ A B C ].T \n");
printMat(ABC_T,3,3);

I= Matmul(ABC,ABC_T,3,3,3); // ABC*ABC.T
printf("\n Identity matrix \n");
printMat(I,3,3);

double u1,u2,u3;
u1=norm(Matcol(I,3,0));// norm of 0th column of ABC
u2=norm(Matcol(I,3,1));// norm of 1th column of ABC
u3=norm(Matcol(I,3,2));// norm of 2nd column of ABC
printf("\nu1= %lf, u2= %lf, u3=%lf\n",u1,u2,u3);

if(round(u1)==round(u2)==round(u3)){
	D=Matadd(A,Matadd(B,C,3,1),3,1); // D= A+B+C
	printf("\n D\n");
	printMat(D,3,1);
        double D_norm,C_norm,B_norm,A_norm;
	D_norm=norm(D);
	C_norm=norm(C);
	B_norm=norm(B);
	A_norm=norm(A);
	double theta1,theta2,theta3;
	theta1=Angle(A,D,A_norm,D_norm);
	printf("\nTheta1=%lf\n",theta1);
	theta2=Angle(B,D,B_norm,D_norm);
	printf("\nTheta2=%lf\n",theta2);
	theta3=Angle(C,D,C_norm,D_norm);
	printf("\nTheta3=%lf\n",theta3);
	printf("\n A+B+c equally inclined to A,B,C\n");
}
else{
	printf("\n A+B+C are not equally inclined to A, B, C\n");
}
}
