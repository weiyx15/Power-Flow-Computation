/*************************FLOW.H****************************/
#ifndef _FLOW_H
#define _FLOW_H
#include <stdio.h>
#define f1(i) (i-1)  
/* 把习惯的一阶矩阵的下标转化为C语言数组下标*/

#define f2(i,j,n) ((i-1)*(n)+j-1)
/* 把习惯的二阶矩阵的下标转化为C语言数组下标*/

FILE *file2=NULL,*file4=NULL,*file6=NULL;

void ybus(int n,int l,int m,float *g,float *b,float *g1,float *b1,float *c1,\
          float *c,float *co,int k,int *s1,int *e1);
void dpqc(float *p,float *q,float *p0,float *q0,float *v,float *v0,int m,\
          int n,float *e,float *f,int k,float *g,float *b,float *dd);
void jmcc(int m,int n,int n0,float *e,float *f,float *g,float *b,float *jm,int k);
void sevc ( float a[], int n0, int k, int n1);
void plsc(int n,int l,int m,float g[],float b[],float e[],float f[],\
           int e1[],int s1[],float g1[],float b1[],float c1[],float c[],\
           float co[],float p1[],float q1[],float p2[],float q2[],float p3[],\
           float q3[],float p[],float q[],float v[],float angle[],int k1);
#endif