
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BINOMIAL	// �Q���������g�p����ێw�肷��

	/* ���[0, 1]�̈�l���� */
#define Rand()		((double)rand()/RAND_MAX)
	/* ���[0, n-1]�̐������� */
#define RandInt(n)	((int)(rand()/(RAND_MAX+1.0)*(n)))

	/* �̂̍\���� */
	typedef struct {
	double *x;	/* genes */
	double f;	/* fitness value */
} IndividualRec, *Individual;

/*
Differential Evolution�̃p�����[�^
*/

/* ���s��(���ϐ��\�𒲂ׂ邽��) */
#define RUN_MAX		50
/* �̐� */
#define Nindividuals	50
/* �ő吢�㐔(�J��Ԃ���) */
#define T_MAX		2000
/* �X�P�[�����O�t�@�N�^�[ */
#define F		0.5
/* ������ */
#define CR		0.5

/*
���̒�`
*/
#define SCHWEFEL221
/* �㉺������̒l */
#define LOWER	-100
#define UPPER	100

/* ����ϐ��x�N�g���̎��� */
#define Nvariables	30

/* �ŏ����C�ő剻�̎w��(���L�͍ŏ���) */
#define better(y1, y2)	(y1<y2)

/* The following is the function of Sum_i (x_i-1)^2 */
void Evaluate(Individual P)
{
	int i;

#ifdef SCHWEFEL221
	// �ő��Βl(Schwefel 2.21)
	P->f = fabs(P->x[0]);
	for (i = 1; i<Nvariables; i++) {
		double v = fabs(P->x[i]);
		if (v>P->f) P->f = v;
	}
#else
	// ���ʊ֐�
	P->f = 0.0;
	for (i = 0; i<Nvariables; i++)
		P->f += (P->x[i] - 1)*(P->x[i] - 1);
#endif
}

/* �̂̏������ƕ]�� */
/* �ŗǌ̂̃C���f�b�N�X��Ԃ� */
int Initialize(Individual P, int n)
{
	int i, j;
	int best;	// �ŗǌ̂̃C���f�b�N�X

	best = 0;
	for (i = 0; i<n; i++) {
		for (j = 0; j<Nvariables; j++)
			P[i].x[j] = LOWER + (UPPER - LOWER)*Rand();
		// �㉺��������Ƀ����_���ɐ���
		Evaluate(&P[i]);	// �̂̕]��
		if (better(P[i].f, P[best].f)) best = i;
	}
	return best;
}

/* �C�ӂ̃f�[�^�^�𓮓I�Ɋ��蓖�Ă�}�N���Ɗ֐� */
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

/* n�̂𓮓I�Ɋ��蓖�Ă�֐� */
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

/* �̂��R�s�[����֐� */
void CopyIndividual(Individual dest, Individual src)
{
	int j;

	for (j = 0; j<Nvariables; j++)
		dest->x[j] = src->x[j];
	dest->f = src->f;
}

/* �̂��o�͂���֐� */
void Print(Individual P)
{
	int j;

	for (j = 0; j<Nvariables; j++)
		printf("%f ", P->x[j]);
	printf(" = %g\n", P->f);
}

/* �ˑR�ψقƌ��� */
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

/* �P���ȍ����i���A���S���Y�� */
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
		srand((run + 1) * 123456789);	// �����V�[�h�̐ݒ�
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