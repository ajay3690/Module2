#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libs/listgen.h"
#include "libs/listfun.h"

int  main()
{
avyuh *A,*B,*C,*D,*E,*F,*v1,*v2,*v3,*vert,*b1,*b11,*v3_a,*v3_b,*c1,*c11,*abc,*abc_t;
avyuh *I;
int m =3,n=1,k=3; //(mxn) matricesavyuh
vert = loadList("Vec.dat",m,k);//load data from dat file
printList(vert);
v1 = Listcol(vert,0);
v2 = Listcol(vert,1);
v3= Listcol(vert,2);

double vn1,vn2,vn3;
vn1=Listnorm(v1);
A=Listscale(v1,1/vn1);
printf("\nA \n");
printList(A);

b1=tr_mul_mat(v2,A,A);
b11=Listsub(v2,b1);
vn2=Listnorm(b11);
B=Listscale(b11,1/vn2);
printf("\nB \n");
printList(B);

v3_a=tr_mul_mat(v3,A,A);//(v3.T @ A )* A
v3_b=tr_mul_mat(v3,B,B);//(v3.T @  B)* B
c11=Listsub(v3,v3_a);
c1=Listsub(c11,v3_b);
vn3=Listnorm(c1);
//printf("\n norm of c1 vn3 = %lf\n\n ",vn3);
C= Listscale(c1,1/vn3);
printf("\n C \n");
printList(C);

abc=VertToList(transposeList(A),transposeList(B),transposeList(C));// [A B C]
printf("\n [ A B C ] \n");
printList(abc);

abc_t=transposeList(abc);
printf("\n tranpose([ A B C ]) \n");
printList(abc_t);

I= Listmul(abc,abc_t);// ABC*ABC.T
printf("\n Identity matrix\n");
printList(I);

double u1, u2, u3;

u1= Listnorm(Listcol(I,0));
u2= Listnorm(Listcol(I,1));
u3= Listnorm(Listcol(I,2));

printf("\nu1= %lf, u2= %lf, u3=%lf\n",u1,u2,u3);

if (round(u1)==round(u2)==round(u3)){
	avyuh *D;
	D=Listadd(A,Listadd(B,C));// D= A+B+C
	printf("\n D \n ");
	printList(D);
        double D_norm,C_norm,B_norm,A_norm;
	D_norm=Listnorm(D);
	C_norm=Listnorm(C);
	B_norm=Listnorm(B);
	A_norm=Listnorm(A);
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
