
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BINOMIAL	// ２項交叉を使用する際指定する

	/* 区間[0, 1]の一様乱数 */
#define Rand()		((double)rand()/RAND_MAX)
	/* 区間[0, n-1]の整数乱数 */
#define RandInt(n)	((int)(rand()/(RAND_MAX+1.0)*(n)))

	/* 個体の構造体 */
	typedef struct {
	double *x;	/* genes */
	double f;	/* fitness value */
} IndividualRec, *Individual;

/*
Differential Evolutionのパラメータ
*/

/* 試行回数(平均性能を調べるため) */
#define RUN_MAX		50
/* 個体数 */
#define Nindividuals	50
/* 最大世代数(繰り返し回数) */
#define T_MAX		2000
/* スケーリングファクター */
#define F		0.5
/* 交叉率 */
#define CR		0.5

/*
問題の定義
*/
#define SCHWEFEL221
/* 上下限制約の値 */
#define LOWER	-100
#define UPPER	100

/* 決定変数ベクトルの次元 */
#define Nvariables	30

/* 最小化，最大化の指定(下記は最小化) */
#define better(y1, y2)	(y1<y2)

/* The following is the function of Sum_i (x_i-1)^2 */
void Evaluate(Individual P)
{
	int i;

#ifdef SCHWEFEL221
	// 最大絶対値(Schwefel 2.21)
	P->f = fabs(P->x[0]);
	for (i = 1; i<Nvariables; i++) {
		double v = fabs(P->x[i]);
		if (v>P->f) P->f = v;
	}
#else
	// 球面関数
	P->f = 0.0;
	for (i = 0; i<Nvariables; i++)
		P->f += (P->x[i] - 1)*(P->x[i] - 1);
#endif
}

/* 個体の初期化と評価 */
/* 最良個体のインデックスを返す */
int Initialize(Individual P, int n)
{
	int i, j;
	int best;	// 最良個体のインデックス

	best = 0;
	for (i = 0; i<n; i++) {
		for (j = 0; j<Nvariables; j++)
			P[i].x[j] = LOWER + (UPPER - LOWER)*Rand();
		// 上下限制約内にランダムに生成
		Evaluate(&P[i]);	// 個体の評価
		if (better(P[i].f, P[best].f)) best = i;
	}
	return best;
}

/* 任意のデータ型を動的に割り当てるマクロと関数 */
#define New(type, n, msg)	(type *)NewCell(sizeof(type), n, msg)

void *NewCell(int size, int n, char *msg)
{
	void *new;

	if ((new = malloc(size*n)) == NULL) {
		fprintf(stderr, "Cannot allocate memory for %d %s\n", n, msg);
		exit(1);
	}
	return new;
}

/* n個体を動的に割り当てる関数 */
Individual NewIndividuals(int n)
{
	int i;
	Individual P;

	P = New(IndividualRec, n, "individuals");
	for (i = 0; i<n; i++) {
		P[i].x = New(double, Nvariables, "x");
	}
	return P;
}

/* 個体をコピーする関数 */
void CopyIndividual(Individual dest, Individual src)
{
	int j;

	for (j = 0; j<Nvariables; j++)
		dest->x[j] = src->x[j];
	dest->f = src->f;
}

/* 個体を出力する関数 */
void Print(Individual P)
{
	int j;

	for (j = 0; j<Nvariables; j++)
		printf("%f ", P->x[j]);
	printf(" = %g\n", P->f);
}

/* 突然変異と交叉 */
void DEoperation(Individual New, Individual Old, int i)
{
	int p1, p2, p3, j, l;

	do {
		p1 = RandInt(Nindividuals);
	} while (p1 == i);
	do {
		p2 = RandInt(Nindividuals);
	} while (p2 == i || p2 == p1);
	do {
		p3 = RandInt(Nindividuals);
	} while (p3 == i || p3 == p1 || p3 == p2);
#ifdef BINOMIAL
	/* DE/rand/1/bin */
	j = RandInt(Nvariables);
	for (l = 0; l<Nvariables; l++) {
		if (l == j || Rand()<CR)
			New[i].x[l] = Old[p1].x[l] + F*(Old[p2].x[l] - Old[p3].x[l]);
		else
			New[i].x[l] = Old[i].x[l];
	}
#else
	/* DE/rand/1/exp */
	CopyIndividual(&New[i], &Old[i]);
	j = RandInt(Nvariables);
	l = 0;
	do {
		New[i].x[j] = Old[p1].x[j] + F*(Old[p2].x[j] - Old[p3].x[j]);
		j = (j + 1) % Nvariables;
		l++;
	} while (l<Nvariables && Rand()<CR);
#endif
}

/* 単純な差分進化アルゴリズム */
void SDE(Individual Pop, Individual New, Individual Best)
{
	int t, i;
	int best, p1, p2;
	double delta;
	Individual Temp;

	best = Initialize(Pop, Nindividuals);
	CopyIndividual(Best, &Pop[best]);

	for (t = 1; t <= T_MAX; t++) {
		for (i = 0; i<Nindividuals; i++) {
			DEoperation(New, Pop, i);
			Evaluate(&New[i]);
			if (better(New[i].f, Pop[i].f)) {
				if (better(New[i].f, Best->f))
					CopyIndividual(Best, &New[i]);
			}
			else
				CopyIndividual(&New[i], &Pop[i]);
		}
		Temp = Pop; Pop = New; New = Temp;
	}
}

int main(void)
{
	int run;
	Individual Pop, New, Best;
	double min, max;

	Pop = NewIndividuals(Nindividuals);
	New = NewIndividuals(Nindividuals);
	Best = NewIndividuals(1);

	double sum = 0, ssum = 0;
	for (run = 0; run<RUN_MAX; run++) {
		srand((run + 1) * 123456789);	// 乱数シードの設定
		SDE(Pop, New, Best);
		printf("%2d: ", run + 1); Print(Best);
		sum += Best->f;
		ssum += Best->f*Best->f;
		if (run == 0 || Best->f<min) min = Best->f;
		if (run == 0 || Best->f>max) max = Best->f;
	}
	double avg = sum / RUN_MAX;
	double std = sqrt(ssum / RUN_MAX - avg*avg);
	printf("Average=%g, Std=%g, Min=%g, Max=%g\n", avg, std, min, max);
	return 0;
}