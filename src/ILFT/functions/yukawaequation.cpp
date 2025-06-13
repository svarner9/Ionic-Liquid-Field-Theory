#include "../ILFT.h"

void ILFT::solveYukawa()
{

	calcChargeDensity();

	double *A, *B, *C, *E, *L, *R, *Q;
	double dr, dr2, coe, p, q, YAdjust; //kappa2, bi;

	A = new double[M];
	B = new double[M];
	C = new double[M];
	E = new double[M];
	L = new double[M];
	R = new double[M - 1];
	Q = new double[M];

	// Constants
	//kappa2 = kappa * kappa;
	dr = gridSize;
	dr2 = gridSize * gridSize;
	coe = 4.0 * Pi * BJ * dr2 * alpha;
	p = 2 + dr2 / (ell * ell) + 2 * dr / ell;
	q = 2 + dr2 / (ell * ell);

	// Matrix entries at boundaries
	A[0] = 0.0;
	B[0] = -p;
	C[0] = 2.0;
	E[0] = rhoZ[0] * coe;
	A[M - 1] = 2.0;
	B[M - 1] = -p;
	C[M - 1] = 0.0;
	E[M - 1] = rhoZ[M - 1] * coe;

	// Fill in matrix interior
	for (int j = 1; j < M - 1; j++)
	{
		A[j] = 1.0;
		B[j] = -q; // - coe;
		C[j] = 1.0;
		E[j] = rhoZ[j] * coe; // - coe*Y[j];
	}

	// Chasing method for tridiagonal matrix
	L[0] = B[0];
	Q[0] = E[0] / L[0];
	for (int j = 1; j < M; j++)
	{
		R[j - 1] = C[j - 1] / L[j - 1];
		L[j] = B[j] - A[j] * R[j - 1];
		Q[j] = (E[j] - A[j] * Q[j - 1]) / L[j];
	}

	// Back substitution
	Y[M - 1] = Q[M - 1];
	for (int j = M - 2; j >= 0; j--)
	{
		Y[j] = Q[j] - R[j] * Y[j + 1];
	}

	// Enforce symmetry
    for (int j = 0; j < M/2; j++)
    {
        Y[M-(j+1)] = -Y[j];
    }

	// YAdjust = -Y[M / 2] + 0.0;
	// for (int j = 0; j < M; j++)
	// {
	// 	Y[j] += YAdjust;
	// }

	delete[] A;
	delete[] B;
	delete[] C;
	delete[] E;
	delete[] L;
	delete[] R;
	delete[] Q;
}


void ILFT::solveYukawa2()
{

    calcChargeDensity();

    double *A, *B, *C, *E, *L, *R, *Q;
    double dr, dr2, coe, p, q, YAdjust; //kappa2, bi;

    A = new double[M];
    B = new double[M];
    C = new double[M];
    E = new double[M];
    L = new double[M];
    R = new double[M - 1];
    Q = new double[M];

    // Constants
    //kappa2 = kappa * kappa;
    dr = gridSize;
    dr2 = gridSize * gridSize;
	coe = 4.0 * Pi * BJ * dr2 * alpha;
    p = 2 + dr2/(ell*ell) + 2*dr/ell;
    q = 2 + dr2/(ell*ell);

    // Matrix entries at boundaries
    A[0] = 0.0;
    B[0] = -p;
    C[0] = 2.0;
    E[0] = rhoZ[0] * coe;
    A[M - 1] = 2.0;
    B[M - 1] = -p;
    C[M - 1] = 0.0;
    E[M - 1] = rhoZ[M-1] * coe;

    // Fill in matrix interior
    for (int j = 1; j < M - 1; j++)
    {
        A[j] = 1.0;
        B[j] = -q;// - coe;
        C[j] = 1.0;
        E[j] = rhoZ[j] * coe;// - coe*Y[j];
    }

    // Chasing method for tridiagonal matrix
    L[0] = B[0];
    Q[0] = E[0] / L[0];
    for (int j = 1; j < M; j++)
    {
        R[j - 1] = C[j - 1] / L[j - 1];
        L[j] = B[j] - A[j] * R[j - 1];
        Q[j] = (E[j] - A[j] * Q[j - 1]) / L[j];
    }

    // Back substitution
    YNew[M - 1] = Q[M - 1];
    for (int j = M - 2; j >= 0; j--)
    {
        YNew[j] = Q[j] - R[j] * YNew[j + 1];
    }

	YAdjust = -YNew[M / 2] + 0.0;
	for (int j = 0; j < M; j++)
	{
		YNew[j] += YAdjust;
	}

	delete[] A;
    delete[] B;
    delete[] C;
    delete[] E;
    delete[] L;
    delete[] R;
    delete[] Q;
}

// void ILFT::solveYukawaf()
// {

//     calcChargeDensityf();

//     double *A, *B, *C, *E, *L, *R, *Q;
//     double dr, dr2, bi, kappa2, coe, p, q;

//     A = new double[M];
//     B = new double[M];
//     C = new double[M];
//     E = new double[M];
//     L = new double[M];
//     R = new double[M - 1];
//     Q = new double[M];

//     // Constants
//     kappa2 = kappa * kappa;
//     dr = gridSize;
//     dr2 = gridSize * gridSize;
//     bi = dr2 * kappa2;
//     coe = 4.0 * Pi * BJ * dr2 * alpha;
//     p = 2 + dr2 / (ell * ell) + 2 * dr / ell;
//     q = 2 + dr2 / (ell * ell);

//     // Matrix entries at boundaries
//     A[0] = 0.0;
//     B[0] = -p;
//     C[0] = 2.0;
//     E[0] = rhoZf[0] * coe;
//     A[M - 1] = 2.0;
//     B[M - 1] = -p;
//     C[M - 1] = 0.0;
//     E[M - 1] = rhoZf[M - 1] * coe;

//     // Fill in matrix interior
//     for (int j = 1; j < M - 1; j++)
//     {
//         A[j] = 1.0;
//         B[j] = -q - bi;
//         C[j] = 1.0;
//         E[j] = rhoZf[j] * coe - bi*Yf[j];
//     }

//     // Chasing method for tridiagonal matrix
//     L[0] = B[0];
//     Q[0] = E[0] / L[0];
//     for (int j = 1; j < M; j++)
//     {
//         R[j - 1] = C[j - 1] / L[j - 1];
//         L[j] = B[j] - A[j] * R[j - 1];
//         Q[j] = (E[j] - A[j] * Q[j - 1]) / L[j];
//     }

//     // Back substitution
//     Yf[M - 1] = Q[M - 1];
//     for (int j = M - 2; j >= 0; j--)
//     {
//         Yf[j] = Q[j] - R[j] * Yf[j + 1];
//     }


//     delete[] A;
//     delete[] B;
//     delete[] C;
//     delete[] E;
//     delete[] L;
//     delete[] R;
//     delete[] Q;
// }