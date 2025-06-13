#include "../ILFT.h"

void ILFT::calcChargeDensity()
{
    for (int j = 0; j < M; j++)
    {
        rhoZ[j] = Z[0] * rho[0][j] + Z[1] * rho[1][j];
    }
}


void ILFT::calculateError()
{
    double temp = 0.0;

    // Set the error from the last iteration
    error_prev = error;

    // Calculate error
    error = 0.0;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < M; j++)
        {
            temp = fabs((rhoNew[i][j] - rho[i][j]) / rhoB[i]);

            if (temp > error)
                error = temp;
        }
    }
}

void ILFT::checkChargeNeutrality()
{
    double out;
    double *rhoZnew;
    rhoZnew = new double[M];

    for (int j = 0; j < M; j++)
    {
        rhoZnew[j] = Z[0] * rhoNew[0][j] + Z[1] * rhoNew[1][j];
    }

    out = simps1(rhoZnew, 0, M-1) * gridSize + surfBC[0] + surfBC[1];
    std::cout << " Net Charge  = " << out << std::endl;
    delete[] rhoZnew;
}

void ILFT::calcNetCharge()
{
    netcharge = trapz(rhoZ,0,M-1) * gridSize + surfCharge[0] + surfCharge[1]; // change to trapz
}

void ILFT::calcSurfaceCharge()
{
	surfCharge[0] = -1*(1/(4*Pi*BJ))*(Psi[1]-Psi[0])/gridSize - gridSize*rhoZ[0]/2;
	surfCharge[1] = 1/(4*Pi*BJ)*(Psi[M-1]-Psi[M-2])/gridSize - gridSize*rhoZ[M-1]/2;
}

int ILFT::getsgn(double x)
{

    int val = 0;

    if (x < 0)
    {
        val = -1;
    }
    else if (x > 0)
    {
        val = 1;
    }
    else if (x == 0)
    {
        val = 0;
    }
    else
    {
        std::cout << "Error in getsgn function" << " \n";
    }

    return val;
}

void ILFT::calcDensityInitial()
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < M; j++)
        {
            rho[i][j] = exp(-Z[i] * (Psi[j] + Y[j])) / ((1 / rhoB[i] - 2)*exp(-1*H[i][j]) + 2 * cosh(Psi[j] + Y[j]));
        }
    }
}

void ILFT::calcChargeDensityf()
{
    for (int j = 0; j < M; j++)
    {
        rhoZf[j] = Z[0] * rhof[0][j] + Z[1] * rhof[1][j];
    }
}

void ILFT::calcDensity(double x)
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < M; j++)
        {
            rhof[i][j] = exp(-Z[i] * (Psi[j] + x + Y[j])) / ((1 / rhoB[i] - 2) * exp(-1 * H[i][j]) + 2 * cosh((Psi[j] + x + Y[j])));
        }
    }
}