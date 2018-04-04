 /*     ���������õı���˵������:                                  *
 *       N:����ڵ�����.         M:�����PQ�ڵ���.                 *
 *       L:�����֧·����.       N0:�ſɱȾ��������.              *
 *       N1:N0+1                 K:��ӡ����.K=1,���ӡ;����,����ӡ.*
 *       K1:�ӳ���PLSC���ж������ѹ����ʽ.K1=1,��Ϊ��������ʽ.����*
 *          Ϊֱ��������ʽ.                                        *
 *       D:�й����޹������������ֵ������ֵ��.                   *
 *       G(I,J):Ybus�ĵ絼Ԫ��(ʵ��).                              *
 *       B(I,J):Ybus�ĵ���Ԫ��(�鲿).                              *
 *       G1(I) :��I֧·�Ĵ����絼.      B1(I):��I֧·�Ĵ�������.   *
 *       C1(I) :��I֧·��pie�ͶԳƽӵص���.                        *
 *       C(I,J):��I�ڵ�J֧·���Գƽӵص���.                        *
 *       CO(I) :��I�ڵ�Ľӵص���.                                 *
 *       S1(I) :��I֧·����ʼ�ڵ��.    E1(I):��I֧·����ֹ�ڵ��. *
 *       P(I)  :��I�ڵ��ע���й�����.  Q(I):��I�ڵ��ע���޹�����.*
 *       P0(I) :��I�ڵ��й��������.    Q0(I):��I�ڵ��޹��������. *
 *       V0(I) :��I�ڵ�(PV�ڵ�)�ĵ�ѹ���(ƽ�����).               *
 *       V(I)  :��I�ڵ�ĵ�ѹ��ֵ.                                 *
 *       E(I)  :��I�ڵ�ĵ�ѹ��ʵ��.    F(I):��I�ڵ�ĵ�ѹ���鲿.  *
 *      JM(I,J):Jacoby����ĵ�I��J��Ԫ��.                          *
 *       A(I,J):�������̵��������,���ǻ�����ĵ�I��J��Ԫ��,����� *
 *              ����A��������һ�д�������Ľ�.                   *
 *       P1(I) :��I֧·��S1(I)�ڵ�ע����й�����.                  *
 *       Q1(I) :��I֧·��S1(I)�ڵ�ע����޹�����.                  *
 *       P2(I) :��I֧·��E1(I)�ڵ�ע����й�����.                  *
 *       Q2(I) :��I֧·��E1(I)�ڵ�ע����޹�����.                  *
 *       P3(I) :��I֧·���й��������.                             *
 *       Q3(I) :��I֧·���޹��������.                             *
 *     ANGLE(I):��I�ڵ��ѹ�ĽǶ�.                                 *
 *     �ڵ���˳��N���ڵ��У�ǰM��ΪPQ�ڵ㣬M+1��N-1���ڵ�Ϊ    *
 *                   PV�ڵ㣬��N��Ϊƽ��ڵ�                       *
*******************************************************************/

// ����2 PQ�ڵ��ת��ΪPQ�ڵ�İ汾
#include "FLOW.H"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <vector>
#include <fstream>
using namespace std;

#define EPSILON 1e-5								// ���������о�

int main()
{
	int N, M, L, K, K1;	// �ڵ�����PQ�ڵ�����֧·��, ��ӡ���أ������꿪��
	float D = 1;					// �й����޹������������ֵ������ֵ��
	float *dd = &D;					// dd: ָ��D��ָ��
	errno_t err;					// �������
	FILE *fid;
	err= fopen_s(&fid,"data_2.txt","r");				// ���������ļ�data_2.txt
	if (err!=0)											// ������룡=0�����ִ���
	{
		exit(1);
	}

	fscanf_s(fid, "%d%d%d%d%d", &N, &M, &L, &K, &K1);
	int N0 = 2*(N-1);								// �ſɱȾ�������
	int N1 = N0+1;									// N0+1
	float *G = (float*) malloc(N*N*sizeof(float));	// �ڵ�絼����
	float *B = (float*) malloc(N*N*sizeof(float));	// �ڵ���ɾ���
	float *G1 = (float*) malloc(L*sizeof(float));	// ֧·�����絼
	float *B1 = (float*) malloc(L*sizeof(float));	// ֧·��������
	float *C1 = (float*) malloc(L*sizeof(float));	// ֧·pie�ͶԳƽӵص���
	float *C = (float*) malloc(N*L*sizeof(float));	// �ڵ�-֧·���Գƽӵص���
	float *C0 = (float*) malloc(N*sizeof(float));	// �ڵ�ӵص���
	int *S1 = (int*) malloc(L*sizeof(int));			// ֧·��ʼ�ڵ��
	int *E1 = (int*) malloc(L*sizeof(int));			// ֧·��ֹ�ڵ��
	float *P = (float*) malloc(N*sizeof(float));	// �ڵ�ע���й�����
	float *Q = (float*) malloc(N*sizeof(float));	// �ڵ�ע���޹�����
	float *P0 = (float*) malloc(N*sizeof(float));	// �ڵ��й��������
	float *Q0 = (float*) malloc(N*sizeof(float));	// �ڵ��޹��������
	float *V0 = (float*) malloc(N*sizeof(float));	// �ڵ��ѹ��ֵƽ�����
	float *V = (float*) malloc(N*sizeof(float));	// �ڵ��ѹ��ֵ
	float *E = (float*) malloc(N*sizeof(float));	// �ڵ��ѹʵ��
	float *F = (float*) malloc(N*sizeof(float));	// �ڵ��ѹ�鲿
	float *JM = (float*) malloc(N0*N0*sizeof(float));// �ſɱȾ���
	float *A = (float*) malloc(N0*N1*sizeof(float));// �������
	float *P1 = (float*) malloc(L*sizeof(float));	// ��I֧·��S1(I)�ڵ�ע����й�����
	float *Q1 = (float*) malloc(L*sizeof(float));	// ��I֧·��S1(I)�ڵ�ע����޹�����
	float *P2 = (float*) malloc(L*sizeof(float));	// ��I֧·��E1(I)�ڵ�ע����й�����
	float *Q2 = (float*) malloc(L*sizeof(float));	// ��I֧·��E1(I)�ڵ�ע����޹�����
	float *P3 = (float*) malloc(L*sizeof(float));	// ��I֧·���й��������
	float *Q3 = (float*) malloc(L*sizeof(float));	// ��I֧·���޹��������
	float *ANGLE = (float*) malloc(N*sizeof(float));// �ڵ��ѹ�Ƕ�
	vector<int> PQ;									// PQ�ڵ����
	vector<int> PV;									// PV�ڵ����
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

	ofstream fout("iteration.txt");					// ����������ļ�iteration.txt

	ybus(N,L,M,G,B,G1,B1,C1,C,C0,K,S1,E1);			// ���ɾ���
	dpqc(P,Q,P0,Q0,V,V0,M,N,E,F,K,G,B,dd,PQ);		// ��deltaP, deltaQ, deltaV2, ����������ļ�
	int cnt = 0;									// ����������
	vector<int>::iterator it;						// PV�ĵ�����
	while ((*dd) > EPSILON)
	{
		cnt++;
		jmcc(M,N,N0,E,F,G,B,JM,K,PQ,PV);					// ���ſɱȾ���
		for (i=1; i<=N0; i++)						// ���������: JdeltaX = -F(X)
		{
			for (j=1; j<=N0; j++)
			{
				A[f2(i,j,N1)] = JM[f2(i,j,N0)];		// ��������ǰn��Ϊ�ſɱȾ���
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
		sevc(A,N0,K,N1);							// �����Է��������ѹ������
		for (i=1; i<=N-1; i++)
		{
			E[f1(i)] += A[f2(2*i-1,N1,N1)];			// E += deltaE
		}
		for (i=1; i<=N-1; i++)
		{
			F[f1(i)] += A[f2(2*i,N1,N1)];			// F += deltaF
		}

		dpqc(P,Q,P0,Q0,V,V0,M,N,E,F,K,G,B,dd,PQ);	// ��deltaP, deltaQ
		plsc(N,L,M,G,B,E,F,E1,S1,G1,B1,C1,C,C0,P1,Q1,P2,Q2,P3,Q3,P,Q,V,ANGLE,K1,PV);// ��P,Q
		for (i=5; i<=7; i++)
		{
			it = find(PV.begin(), PV.end(), i);
			if (it!=PV.end())		// Q�����޶�
			{
				if (Q[i-1] < -0.1)
				{
					Q[i-1] = -0.1;
					PV.erase(it);
					PQ.push_back(i);
					printf("\n PV�ڵ�%dת��ΪPQ�ڵ�\n",i);
				}
				else if (Q[i-1] > 0.1)
				{
					Q[i-1] = 0.1;
					PV.erase(it);
					PQ.push_back(i);
					printf("\n PV�ڵ�%dת��ΪPQ�ڵ�\n",i);
				}
			}
		}
		printf("\n����������%d, ��%f\n",cnt,(*dd));
		fout << cnt << "\t" << (*dd) << endl;
		system("pause");
	}
	printf("�ܵĵ���������%d",cnt);
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