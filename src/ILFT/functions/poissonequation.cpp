#include "../ILFT.h"

void ILFT::solvePoisson()
{

    // Calculate the charge density
    calcChargeDensity();

    // Decide where to send calculation
    if (BC_type == 0)
    {
        solvePoissonPotential();
    }
    else
    {
        solvePoissonSurfaceCharge();
    }
}

void ILFT::solvePoissonPotential()
{
    // Routine solves Poisson equation when two surfaces have fixed potential

    double *Int1;
    double z, Int2, Int3;
    double psi0, psiL, temp, tempz, coe;

    // Allocate
    Int1 = new double[M];

    // Coefficient
    coe = -4.0 * Pi * BJ;

    psi0 = surfBC[0];
    psiL = surfBC[1];


    // Make array that is z'*rhoZ
    for (int i = 0; i < M; i++)
    {
        z = i * gridSize;
        Int1[i] = z * rhoZ[i];
    }

    // Calculate integral of rhoZ from 0 to L (lower to upper)
    Int2 = trapz(rhoZ, 0, M - 1) * gridSize;
    //Int2 = simps1(rhoZ, 0, M - 1) * gridSize;

    // Calculate integral of z'*rhoZ from 0 to L (lower to upper)
    Int3 = trapz(Int1, 0, M - 1) * gridSize;
    //Int3 = simps1(Int1, 0, M - 1) * gridSize;

    // Solve potential
    for (int i = 0; i < M; i++)
    {
        // Current spatial position
        z = i * gridSize;

        // Calculate integral of rhoZ from 0 to z
        temp = trapz(rhoZ, 0, i) * gridSize;
        //temp = simps1(rhoZ, 0, i) * gridSize;

        // Calculate integral of z'*rhoZ from 0 to z
        tempz = trapz(Int1, 0, i) * gridSize;
        //tempz = simps1(Int1, 0, i) * gridSize;

        // Solve for potential
        Psi[i] = psi0 + (psiL - psi0) * z / separation + coe * z * temp - coe * tempz - coe * z * Int2 + coe * z * Int3 / separation;
    }

    // // Enforce symmetry
    // for (int i = 0; i < (M-1)/2; i++)
    // {
    //     Psi[M-1-i] = -1*Psi[i];
    // }

	// Psi[(M-1)/2] = 0.0;

    delete[] Int1;
}

void ILFT::solvePoissonSurfaceCharge()
{
	double *A, *B, *C, *E, *L, *R, *Q;
	double dr2, coe, PsiAdjust, bi;//, kappa2;

	A = new double[M];
	B = new double[M];
	C = new double[M];
	E = new double[M];
	L = new double[M];
	R = new double[M - 1];
	Q = new double[M];

	// Constants
	// kappa2 = kappa * kappa;
	dr2 = gridSize * gridSize;
	bi = 4.0 * Pi * BJ * dr2;
	coe = 4.0 * Pi * BJ * dr2;

	// Matrix entries at boundaries
	A[0] = 0.0;
	B[0] = -1.0 - bi;
	C[0] = 1.0;
	E[0] = -4.0 * Pi * BJ * surfBC[0] * gridSize - rhoZ[0] * coe / 2 - bi * Psi[0];
	A[M - 1] = 1.0;
	B[M - 1] = -1.0 - bi;
	C[M - 1] = 0.0;
	E[M - 1] = -4.0 * Pi * BJ * surfBC[1] * gridSize - rhoZ[M - 1] * coe / 2 - bi * Psi[M - 1];

	// Fill in matrix interior
	for (int j = 1; j < M - 1; j++)
	{
		A[j] = 1.0;
		B[j] = -2.0 - bi;
		C[j] = 1.0;
		E[j] = -1.0 * rhoZ[j] * coe - bi * Psi[j];
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
	Psi[M - 1] = Q[M - 1];
	for (int j = M - 2; j >= 0; j--)
	{
		Psi[j] = Q[j] - R[j] * Psi[j + 1];
	}

	// Offset Psi
	PsiAdjust = -Psi[M / 2] + xPsi;
	for (int j = 0; j < M; j++)
	{
		Psi[j] += PsiAdjust;
	}

    // Enforce symmetry
    for (int j = 0; j < M/2; j++)
    {
        Psi[M-(j+1)] = -Psi[j];
    }



	delete[] A;
	delete[] B;
	delete[] C;
	delete[] E;
	delete[] L;
	delete[] R;
	delete[] Q;
}

void ILFT::solvePoissonSurfaceCharge2()
{
    double *A, *B, *C, *E, *L, *R, *Q;
    double dr2, coe, PsiAdjust, bi;//, kappa2;

    A = new double[M];
    B = new double[M];
    C = new double[M];
    E = new double[M];
    L = new double[M];
    R = new double[M - 1];
    Q = new double[M];

    // Constants
	// kappa2 = kappa * kappa;
    dr2 = gridSize * gridSize;
	bi = 4.0 * Pi * BJ * dr2;
	coe = 4.0 * Pi * BJ * dr2;

    // Matrix entries at boundaries
    A[0] = 0.0;
    B[0] = -1.0 - bi;
    C[0] = 1.0;
    E[0] = -4.0 * Pi * BJ * surfBC[0] * gridSize - rhoZ[0] * coe / 2 - bi * Psi[0];
    A[M - 1] = 1.0;
    B[M - 1] = -1.0 - bi;
    C[M - 1] = 0.0;
    E[M - 1] = -4.0 * Pi * BJ * surfBC[1] * gridSize - rhoZ[M - 1] * coe / 2 - bi * Psi[M-1];

    // Fill in matrix interior
    for (int j = 1; j < M - 1; j++)
    {
        A[j] = 1.0;
        B[j] = -2.0 - bi;
        C[j] = 1.0;
        E[j] = -1.0 * rhoZ[j] * coe - bi * Psi[j];
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
    PsiNew[M - 1] = Q[M - 1];
    for (int j = M - 2; j >= 0; j--)
    {
        PsiNew[j] = Q[j] - R[j] * PsiNew[j + 1];
    }

    // Offset Psi
    PsiAdjust = -PsiNew[M / 2] + xPsi;
    for (int j = 0; j < M; j++)
    {
        PsiNew[j] += PsiAdjust;
    }

    // Enforce symmetry
    for (int j = 0; j < M/2; j++)
    {
        PsiNew[M-(j+1)] = PsiNew[j];
    }
    
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] E;
    delete[] L;
    delete[] R;
    delete[] Q;
}

/**************************************************/
/* Code adapted from code written by Chris Balzer */
/**************************************************/
void ILFT::chargeNeutral()
{
    /* Root finding for electrostatic potential using Steffensen's method
     Trying to solve...
        int_0^L rhoZ + sigma0 + sigmaL = 0
    */

    int iterin;
    int maxIterations = 1E6;
    double fErrorTol = 1E-9;
    double sigma0, sigmaL;
    double x, f, fh, g;

    double rho1, rho2, rho3, rho4;
    double x1, x2;
    double f1, f2, fh1;

    double *intRho;
    intRho = new double[2];

    // More expressive variable names
    sigma0 = surfBC[0];
    sigmaL = surfBC[1];

    // Initial guess
    x = 0.0;
    
    calcDensity(x);

    // Get contribution to total charge from each component
    for (int i = 0; i < 2; i++)
    {
        intRho[i] = simps1(rhof[i], 0, M-1) * gridSize;
    }

    // Calculate value of root function
    f = sigma0 + sigmaL;
    for (int i = 0; i < 2; i++)
    {
        f += Z[i] * intRho[i];
    }

    iterin = 0;
    while (fabs(f) > fErrorTol && iterin < maxIterations)
    {
        rho1 = rhof[1][0];
        x1 = x;
        f1 = f;

        calcDensity(x+f);

        rho2 = rhof[1][0];

        for (int i = 0; i < 2; i++)
        {
            intRho[i] = trapz(rhof[i], 0, M-1) * gridSize;
        }
        // Calculate one step ahead
        fh = sigma0 + sigmaL;
        for (int i = 0; i < 2; i++)
        {
            fh += Z[i] * intRho[i];
        }

        fh1 = fh;

        // Calculate g
        g = fh / f - 1.0;

        if (g == 0)
        {
            std::cout << "g = " << g << std::endl;
            break;
        }

        // Update the guess for x
        x -= f / g;

        x2 = x;

        calcDensity(x);

        rho3 = rhof[1][0];

        for (int i = 0; i < 2; i++)
        {
            intRho[i] = trapz(rhof[i], 0, M-1) * gridSize;
        }
        // Calculate value of root function with updated guess
        f = sigma0 + sigmaL;
        for (int i = 0; i < 2; i++)
        {
            f += Z[i] * intRho[i];
        }

        f2 = f;

        if (std::isnan(f))
        {
            calcDensity(x1);
            writeRhof();
        }

        // Increase iteration count
        iterin++;

    }

    // Update outer value
    xPsi += x;

    calcDensity(x);

    rho4 = rhof[1][0];

    if (std::isnan(rhof[0][0]-rhof[1][0]) || std::isnan(f))
    {
        std::cout << "Densty isnan!!!!" << std::endl;
        std::cout << "f =        " << f << std::endl;
        std::cout << "f1 =       " << f1 << std::endl;
        std::cout << "f2 =       " << f2 << std::endl;
        std::cout << "fh =       " << fh << std::endl;
        std::cout << "fh1 =      " << fh1 << std::endl;
        std::cout << "g =        " << g << std::endl;
        std::cout << "x =        " << x << std::endl;
        std::cout << "x1 =       " << x1 << std::endl;
        std::cout << "x2 =       " << x2 << std::endl;
        std::cout << "xPsi =     " << xPsi << std::endl;
        std::cout << "initer =   " << iterin << std::endl;
        std::cout << "outiter =  " << iter << std::endl;
        std::cout << "rho1 =     " << rho1 << std::endl;
        std::cout << "rho2 =     " << rho2 << std::endl;
        std::cout << "rho3 =     " << rho3 << std::endl;
        std::cout << "rho4 =     " << rho4 << std::endl;
        writeLog("chargeNeutral() - isnan triggered!");
        finalizeIO();
        exit(0);
    }

    // Update density to ensure charge neutrality
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < M; j++)
        {
            rhoNew[i][j] = rhof[i][j];
        }
    }

    delete[] intRho;

    // Update the potential
    if (fabs(f) > fErrorTol)
    {
        std::cout << "Code failed. See log file." << std::endl;
        writeLog("chargeNeutral() - could not find x!");
        finalizeIO();
        exit(0);
    }
}