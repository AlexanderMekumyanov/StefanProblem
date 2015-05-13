#include <stdio.h>
#include <conio.h>
#include <math.h>
#include "stefanparams.h"
#include "progonka.h"

void   FillMatrix();
void   FillMatrixNew();
double CroEff(int i);
void   CalcSchemaParams();
void   WriteArrayToFile(char *fileName, double *arr, int size);

void CalcTwoPhaseStefanProblem()
{	
	Progonka* progonka = Progonka::GetInstance();

	x[0]  = 0;
	h[1]  = L / (size-1);
	x[1]  = x[0] + h[1];
	hs[1] = h[1]/2;
	for(int i=2;i<100;i++)
	{
		h[i]  = L / (size-1);
		x[i]  = x[i-1] + h[i];
		hs[i] = (h[i] + h[i+1]) / 2;
	}
	////initial conditions
	for(int i = 0; i < size; i++)
	{
		T[i] = T0;
	}

	CalcSchemaParams();
	//starting calculation
	tau = tk / num;
	double steps[num];
	currTime = tau;
	int k = 0;
	while (currTime < tk)
	{
		double err = 1;
		int    n   = 0;
		//save previous layer
		
		while(err > eps)
		{
			for(int i = 0; i < size; i++)
			{
				oldT[i] = T[i];
			}

			FillMatrixNew();

			progonka->CalcAlphaBeta(As,Bs,Cs,Ds,alphaT,betaT,size,0,2);
			T[size-1] = (As[size-1] * betaT[size-2] + Ds[size-1]) / (Cs[size-1] - As[size-1] * alphaT[size-2]);
			progonka->SolveWithProgonka(T,alphaT,betaT,size,0,2,err);
			
			//CalcSchemaParams();
			//calculate an error
			n++;
		}
		if(currTime == tau)
		{
			WriteArrayToFile("1Temp.txt",T,size);
		}
		currTime += tau;
		steps[k]  = n;
		k++;
	}
	WriteArrayToFile("Temp.txt", T, size);
	WriteArrayToFile("inerations.txt",steps,num);
	printf("Calculations successfully finished!");
	getch();
}

void FillMatrix()
{
	//As[0] = 0; Bs[0] = 1; Cs[0] = 0; Ds[0] = Tc;
	As[size-1] = 0; 
	Bs[size-1] = CroEff(size-1) * hs[size-1]                / tau + lamda[size-1] / h[size-1];
	Cs[size-1] = lamda[size-1]                              / h[size-1]; 
	Ds[size-1] = CroEff(size-1) * hs[size-1] * oldT[size-1] / tau;

	for(int i=1;i<size-1;i++)
	{
		As[i] = lamda[i+1]                  / h[i+1];
		Bs[i] = lamda[i+1]                  / h[i+1] + lamda[i]/h[i] + CroEff(i) * hs[i] / tau;
		Cs[i] = lamda[i]                    / h[i];
		Ds[i] = CroEff(i) * hs[i] * oldT[i] / tau;
	}
};
void FillMatrixNew()
{
	As[0] = 0; 
	Bs[0] = lamda[1] / h[1]; 
	Cs[0] = Bs[0] + CroEff(0) * hs[0] / tau+alpha; 
	Ds[0] = CroEff(0) * hs[0] / tau * oldT[0] + alpha * Tc;
	for(int i = 1;i < size - 1;i++)
	{
		As[i] = Bs[i - 1];//lamda[i] / h[i];
		Bs[i] = lamda[i]                          / h[i];
		Cs[i] = As[i] + Bs[i] + CroEff(i) * hs[i] / tau;
		Ds[i] = oldT[i] * CroEff(i)*hs[i]         / tau;
	}
	As[size-1] = Bs[size-2];
	Bs[size-1] = 0;
	Cs[size-1] = As[size-1] + CroEff(size-1) * hs[size-1]   / tau;
	Ds[size-1] = CroEff(size-1) * hs[size-1] * oldT[size-1] / tau;
}
double CroEff(int i)
{
	switch(i)
	{
		case 0:      return m*beta[0]*Crow + (1-m)*Crosc + m*beta[0]*Croi + dP[0];
		case size-1: return m*beta[size-1]*Crow + (1-m)*Crosc+m*beta[size-1]*Croi + dM[size-1];
		default:     return m*beta[i]*Crow + (1-m)*Crosc+m*beta[i]*Croi + dP[i] - dM[i];
	}
}

void CalcSchemaParams()
{
	if(T[0] > Tp)
	{
		beta[0] = 1;
	}
	else
	{
		beta[0] = 0;
	}

	lamda[0] = m * beta[0] * lamdaw + (1-m) * lamdasc + m * beta[0] * lamdai;

	for(int i = 1;i < size; i++)
	{
		if(T[i] > Tp)
		{
			beta[i] = 1;
		}
		else
		{
			beta[i] = 0;
		}
		if(T[i]!=T[i-1])
		{
			teta[i] = (T[i] - Tp)/(T[i] - T[i-1]);
		}
		else
		{
			teta[i] = 0;
		}
		lamda[i] = m * beta[i] * lamdaw + (1-m) * lamdasc + m * beta[i] * lamdai;
	}
	//calc d+ and d-
	if(T[size-1] != T[size-2])
	{
		dM[size-1] = (beta[size-1] - beta[size-2]) *q * row * (teta[size-1]) / (T[size-1]-T[size-2]);
	}
	else 
	{
		dM[size-1] = 0;
	}
	if(T[1] != T[0])
	{
		dP[0] = (beta[1] - beta[0])*q*row*(1-teta[1]) / (T[0]-T[1]);
	}
	else 
	{
		dP[0] = 0;
	}
	for(int i = 1; i < size-1; i++)
	{
		if(T[i] != T[i-1])
		{
			dM[i] = (beta[i] - beta[i-1])*q*row*(teta[i]) / (T[i] - T[i-1]);
		}
		else
		{
			dM[i] = 0;
		}
		if(T[i-1] != T[i])
		{
			dP[i] = (beta[i] - beta[i-1])*q*row*(teta[i]) / (T[i] - T[i-1]);
		}
		else
		{
			dP[i] = 0;
		}
	}
};

void WriteArrayToFile(char *fileName, double *arr, int size)
{
	FILE *f = fopen(fileName, "w");
	for(int i = 0; i < size; i++)
	{
		fprintf(f, "%f\n", arr[i]);
	}
	fclose(f);
};