 /*     程序中所用的变量说明如下:                                  *
 *       N:网络节点总数.         M:网络的PQ节点数.                 *
 *       L:网络的支路总数.       N0:雅可比矩阵的行数.              *
 *       N1:N0+1                 K:打印开关.K=1,则打印;否则,不打印.*
 *       K1:子程序PLSC中判断输入电压的形式.K1=1,则为极座标形式.否则*
 *          为直角坐标形式.                                        *
 *       D:有功及无功功率误差的最大值（绝对值）.                   *
 *       G(I,J):Ybus的电导元素(实部).                              *
 *       B(I,J):Ybus的电纳元素(虚部).                              *
 *       G1(I) :第I支路的串联电导.      B1(I):第I支路的串联电纳.   *
 *       C1(I) :第I支路的pie型对称接地电纳.                        *
 *       C(I,J):第I节点J支路不对称接地电纳.                        *
 *       CO(I) :第I节点的接地电纳.                                 *
 *       S1(I) :第I支路的起始节点号.    E1(I):第I支路的终止节点号. *
 *       P(I)  :第I节点的注入有功功率.  Q(I):第I节点的注入无功功率.*
 *       P0(I) :第I节点有功功率误差.    Q0(I):第I节点无功功率误差. *
 *       V0(I) :第I节点(PV节点)的电压误差(平方误差).               *
 *       V(I)  :第I节点的电压幅值.                                 *
 *       E(I)  :第I节点的电压的实部.    F(I):第I节点的电压的虚部.  *
 *      JM(I,J):Jacoby矩阵的第I行J列元素.                          *
 *       A(I,J):修正方程的增广矩阵,三角化矩阵的第I行J列元素,运算结 *
 *              束后A矩阵的最后一列存放修正的解.                   *
 *       P1(I) :第I支路由S1(I)节点注入的有功功率.                  *
 *       Q1(I) :第I支路由S1(I)节点注入的无功功率.                  *
 *       P2(I) :第I支路由E1(I)节点注入的有功功率.                  *
 *       Q2(I) :第I支路由E1(I)节点注入的无功功率.                  *
 *       P3(I) :第I支路的有功功率损耗.                             *
 *       Q3(I) :第I支路的无功功率损耗.                             *
 *     ANGLE(I):第I节点电压的角度.                                 *
 *     节点编号顺序：N个节点中，前M个为PQ节点，M+1至N-1个节点为    *
 *                   PV节点，第N个为平衡节点                       *
*******************************************************************/

// 算例2 PQ节点可转化为PQ节点的版本
#include "FLOW.H"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <vector>
#include <fstream>
using namespace std;

#define EPSILON 1e-5								// 迭代收敛判据

int main()
{
	int N, M, L, K, K1;	// 节点数，PQ节点数，支路数, 打印开关，极坐标开关
	float D = 1;					// 有功及无功功率误差的最大值（绝对值）
	float *dd = &D;					// dd: 指向D的指针
	errno_t err;					// 错误代码
	FILE *fid;
	err= fopen_s(&fid,"data_2.txt","r");				// 输入数据文件data_2.txt
	if (err!=0)											// 错误代码！=0，出现错误
	{
		exit(1);
	}

	fscanf_s(fid, "%d%d%d%d%d", &N, &M, &L, &K, &K1);
	int N0 = 2*(N-1);								// 雅可比矩阵行数
	int N1 = N0+1;									// N0+1
	float *G = (float*) malloc(N*N*sizeof(float));	// 节点电导矩阵
	float *B = (float*) malloc(N*N*sizeof(float));	// 节点电纳矩阵
	float *G1 = (float*) malloc(L*sizeof(float));	// 支路串联电导
	float *B1 = (float*) malloc(L*sizeof(float));	// 支路串联电纳
	float *C1 = (float*) malloc(L*sizeof(float));	// 支路pie型对称接地电纳
	float *C = (float*) malloc(N*L*sizeof(float));	// 节点-支路不对称接地电纳
	float *C0 = (float*) malloc(N*sizeof(float));	// 节点接地电纳
	int *S1 = (int*) malloc(L*sizeof(int));			// 支路起始节点号
	int *E1 = (int*) malloc(L*sizeof(int));			// 支路终止节点号
	float *P = (float*) malloc(N*sizeof(float));	// 节点注入有功功率
	float *Q = (float*) malloc(N*sizeof(float));	// 节点注入无功功率
	float *P0 = (float*) malloc(N*sizeof(float));	// 节点有功功率误差
	float *Q0 = (float*) malloc(N*sizeof(float));	// 节点无功功率误差
	float *V0 = (float*) malloc(N*sizeof(float));	// 节点电压幅值平方误差
	float *V = (float*) malloc(N*sizeof(float));	// 节点电压幅值
	float *E = (float*) malloc(N*sizeof(float));	// 节点电压实部
	float *F = (float*) malloc(N*sizeof(float));	// 节点电压虚部
	float *JM = (float*) malloc(N0*N0*sizeof(float));// 雅可比矩阵
	float *A = (float*) malloc(N0*N1*sizeof(float));// 增广矩阵
	float *P1 = (float*) malloc(L*sizeof(float));	// 第I支路由S1(I)节点注入的有功功率
	float *Q1 = (float*) malloc(L*sizeof(float));	// 第I支路由S1(I)节点注入的无功功率
	float *P2 = (float*) malloc(L*sizeof(float));	// 第I支路由E1(I)节点注入的有功功率
	float *Q2 = (float*) malloc(L*sizeof(float));	// 第I支路由E1(I)节点注入的无功功率
	float *P3 = (float*) malloc(L*sizeof(float));	// 第I支路的有功功率损耗
	float *Q3 = (float*) malloc(L*sizeof(float));	// 第I支路的无功功率损耗
	float *ANGLE = (float*) malloc(N*sizeof(float));// 节点电压角度
	vector<int> PQ;									// PQ节点序号
	vector<int> PV;									// PV节点序号
	int i = 0, j = 0;
	for (i=1; i<=M; i++)
	{
		PQ.push_back(i);
	}
	for (; i<= N-1; i++)
	{
		PV.push_back(i);
	}



	memset(G, 0, N*N*sizeof(float));
	memset(B, 0, N*N*sizeof(float));
	for (i=0; i<L; i++)
	{
		fscanf_s(fid, "%f", &G1[i]);
	}
	for (i=0; i<L; i++)
	{
		fscanf_s(fid, "%f", &B1[i]);
	}
	for (i=0; i<L; i++)
	{
		fscanf_s(fid, "%f", &C1[i]);
	}
	for (i=0; i<N*L; i++)
	{
		fscanf_s(fid, "%f", &C[i]);
	}
	for (i=0; i<N; i++)
	{
		fscanf_s(fid, "%f", &C0[i]);
	}
	for (i=0; i<L; i++)
	{
		fscanf_s(fid, "%d", &S1[i]);
	}
	for (i=0; i<L; i++)
	{
		fscanf_s(fid, "%d", &E1[i]);
	}
	for (i=0; i<N; i++)
	{
		fscanf_s(fid, "%f", &P[i]);
	}
	for (i=0; i<N; i++)
	{
		fscanf_s(fid, "%f", &Q[i]);
	}
	memset(P0, 0, N*sizeof(int));
	memset(Q0, 0, N*sizeof(int));
	memset(V0, 0, N*sizeof(int));
	for (i=0; i<N; i++)
	{
		fscanf_s(fid, "%f", &V[i]);
	}
	for (i=0; i<N; i++)
	{
		fscanf_s(fid, "%f", &E[i]);
	}
	for (i=0; i<N; i++)
	{
		fscanf_s(fid, "%f", &F[i]);
	}
	memset(JM, 0, N0*N0*sizeof(float));
	memset(A, 0, N0*N1*sizeof(float));
	memset(P1, 0, L*sizeof(float));
	memset(Q1, 0, L*sizeof(float));
	memset(P2, 0, L*sizeof(float));
	memset(Q2, 0, L*sizeof(float));
	memset(P3, 0, L*sizeof(float));
	memset(Q3, 0, L*sizeof(float));
	memset(ANGLE, 0, N*sizeof(float));

	ofstream fout("iteration.txt");					// 输出迭代误差到文件iteration.txt

	ybus(N,L,M,G,B,G1,B1,C1,C,C0,K,S1,E1);			// 求导纳矩阵
	dpqc(P,Q,P0,Q0,V,V0,M,N,E,F,K,G,B,dd,PQ);		// 求deltaP, deltaQ, deltaV2, 输出迭代误差到文件
	int cnt = 0;									// 迭代计数器
	vector<int>::iterator it;						// PV的迭代器
	while ((*dd) > EPSILON)
	{
		cnt++;
		jmcc(M,N,N0,E,F,G,B,JM,K,PQ,PV);					// 求雅可比矩阵
		for (i=1; i<=N0; i++)						// 求增广矩阵: JdeltaX = -F(X)
		{
			for (j=1; j<=N0; j++)
			{
				A[f2(i,j,N1)] = JM[f2(i,j,N0)];		// 增广矩阵的前n列为雅可比矩阵
			}
			if (i%2==0)
			{
				A[f2(i,N1,N1)] = -P0[i/2-1];		// -deltaQ
			}
			else if ((find(PQ.begin(),PQ.end(),i/2+1))!=PQ.end())
			{
				A[f2(i,N1,N1)] = -Q0[i/2];			// -deltaP
			}
			else
			{
				A[f2(i,N1,N1)] = -V0[i/2];			// -deltaV2
			}
		}
		sevc(A,N0,K,N1);							// 解线性方程组求电压修正量
		for (i=1; i<=N-1; i++)
		{
			E[f1(i)] += A[f2(2*i-1,N1,N1)];			// E += deltaE
		}
		for (i=1; i<=N-1; i++)
		{
			F[f1(i)] += A[f2(2*i,N1,N1)];			// F += deltaF
		}

		dpqc(P,Q,P0,Q0,V,V0,M,N,E,F,K,G,B,dd,PQ);	// 求deltaP, deltaQ
		plsc(N,L,M,G,B,E,F,E1,S1,G1,B1,C1,C,C0,P1,Q1,P2,Q2,P3,Q3,P,Q,V,ANGLE,K1,PV);// 求P,Q
		for (i=5; i<=7; i++)
		{
			it = find(PV.begin(), PV.end(), i);
			if (it!=PV.end())		// Q超出限定
			{
				if (Q[i-1] < -0.1)
				{
					Q[i-1] = -0.1;
					PV.erase(it);
					PQ.push_back(i);
					printf("\n PV节点%d转化为PQ节点\n",i);
				}
				else if (Q[i-1] > 0.1)
				{
					Q[i-1] = 0.1;
					PV.erase(it);
					PQ.push_back(i);
					printf("\n PV节点%d转化为PQ节点\n",i);
				}
			}
		}
		printf("\n迭代次数：%d, 误差：%f\n",cnt,(*dd));
		fout << cnt << "\t" << (*dd) << endl;
		system("pause");
	}
	printf("总的迭代次数：%d",cnt);
	plsc(N,L,M,G,B,E,F,E1,S1,G1,B1,C1,C,C0,P1,Q1,P2,Q2,P3,Q3,P,Q,V,ANGLE,K1,PV);

	fout.close();
	free(G);
	free(B);
	free(G1);
	free(B1);
	free(C1);
	free(C);
	free(C0);
	free(S1);
	free(E1);
	free(P);
	free(Q);
	free(P0);
	free(Q0);
	free(V0);
	free(V);
	free(E);
	free(F);
	free(JM);
	free(A);
	free(P1);
	free(Q1);
	free(P2);
	free(Q2);
	free(P3);
	free(Q3);
	free(ANGLE);
	return 0;
}