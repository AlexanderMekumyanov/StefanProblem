#include "progonka.h"
#include <math.h>
#include <stdlib.h>

Progonka* Progonka::instance = NULL;

Progonka::Progonka()
{
}

Progonka::~Progonka()
{
}

Progonka* Progonka::GetInstance()
{
	if (instance != NULL)
		instance = new Progonka();
	return instance;
}

void Progonka::CalcAlphaBeta(double *As, double *Bs, double *Cs, double *Ds, double *alpha, double *beta, int size, int front, int sweepType)
{
	if(sweepType == 0)
	{
		alpha[0] = (Bs[0] - As[0])/Bs[0];
		beta[0]  = Ds[0] / Bs[0];
		for(int i=1;i<front;i++)
		{
			alpha[i] = (Bs[i]- Cs[i] * (1-alpha[i-1]) - As[i]) / (Bs[i]- Cs[i]*(1-alpha[i-1]));
			beta[i]  = (Ds[i]+Cs[i]  * beta[i-1])              / (Bs[i]- Cs[i]*(1-alpha[i-1]));
		}
	}
	else
	{
		if (sweepType == 1)
		{
			alpha[size - 1] = (Bs[size - 1] - Cs[size - 1]) / Bs[size - 1];
			beta[size - 1] = Ds[size - 1] / Bs[size - 1];
			for (int i = size - 2; i > front; i--)
			{
				alpha[i] = (Bs[i] - As[i] * (1 - alpha[i + 1]) - Cs[i]) / (Bs[i] - As[i] * (1 - alpha[i + 1]));
				beta[i] = (Ds[i] + As[i] * beta[i + 1]) / (Bs[i] - As[i] * (1 - alpha[i + 1]));
			}
		}
		else
		{
			alpha[0] = Bs[0] / Cs[0];
			beta[0] = Ds[0] / Cs[0];
			double a, b, c, d, tempAlpha, tempBeta;
			double temp;
			for (int i = 1; i < size - 1; i++)
			{
				if (i == 99)
				{
					a = As[i];
					b = Bs[i];
					c = Cs[i];
					d = Ds[i];
					tempAlpha = alpha[i - 1];
					tempBeta = beta[i - 1];
				}
				temp = Cs[i] - As[i] * alpha[i - 1];
				alpha[i] = Bs[i] / temp;
				beta[i] = (As[i] * beta[i - 1] + Ds[i]) / temp;
			}
		}
	}
}

int Progonka::SolveWithProgonka(double *x, double *alpha, double *beta, int size, int front, int sweepType, double &err)
{
	if(x == NULL || alpha == NULL || beta == NULL || front < 0 || front > size)
	{
		return -1;
	}
	err=0;
	double temp;
//	CalcAlphaBeta(As,Bs,Cs,Ds,alpha,beta,size);
	if(sweepType == 0)
	{
		for(int i=front-1;i>=0;i--)
		{
			temp = x[i];
			x[i] = (1-alpha[i])*x[i+1]+beta[i];
			if(abs(temp - x[i]) > err) err = abs(temp - x[i]);
		}
	}
	else if(sweepType == 1)
	{
		for(int i=front+1;i<size;i++)
		{
			temp = x[i];
			x[i] = (1-alpha[i])*x[i-1]+beta[i];
			if(abs(temp - x[i]) > err) err = abs(temp - x[i]);
		}
	}
	else 
	{
//		x[size-1] = (As[size-1] * beta[size-2] + Ds[size-1]) / (Cs[size-1] - As[size-1] * alpha[size-2]);
		for(int i=size-2;i>=0;i--)
		{
			temp = x[i];
			x[i] = (/*1-*/alpha[i]) * x[i+1] + beta[i];
			if(abs(temp - x[i]) > err) 
				err = abs(temp-x[i]);
		}
	}
	return 0;
}