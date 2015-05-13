// стр 118
class Progonka
{
public:
	Progonka();
	~Progonka();
	
	void   CalcAlphaBeta(double *As, double *Bs, double *Cs, double *Ds, double *alpha, double *beta, int size, int front, int sweepType);
	int    SolveWithProgonka(double *x, double *alpha, double *beta, int size, int front, int sweepType, double &err);

	static Progonka* GetInstance();	
private:
	static Progonka *instance;
};