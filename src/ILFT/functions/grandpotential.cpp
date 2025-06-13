#include "../ILFT.h"

void ILFT::calcGrandPotential()
{
    double Int1, Int2, Int3, Int4, Int5=0;
    double *temp;
    double dPdz, coe, mu_b;

    temp = new double[M];
    coe = 4*Pi*BJ;

    // Entropy Integrand
    for (int j = 0; j < M; j++)
    {
        if (rhoB[0]<0.4999999)
        {
            temp[j] = rho[0][j] * log(rho[0][j]) + rho[1][j] * log(rho[1][j]) + (1 - rho[0][j] - rho[1][j]) * log(1 - rho[0][j] - rho[1][j]);
        }
        else
        {
            temp[j] = rho[0][j] * log(rho[0][j]) + rho[1][j] * log(rho[1][j]);
        }
        
    }

    // Calculate entropy contribution
    Int1 = trapz(temp, 0, M - 1) * gridSize;

    // Electrostatic Integrand
    for (int j = 0; j < M; j++)
    {
        dPdz = (j==0)*(Psi[j+1]-Psi[j])/gridSize + (j==(M-1))*(Psi[j]-Psi[j-1])/gridSize + (j!=0)*(j!=(M-1))*(Psi[j+1]-Psi[j-1])/2/gridSize;
        temp[j] = rhoZ[j]*Psi[j] - 0.5/coe*dPdz*dPdz;
    }

    Int2 = trapz(temp, 0, M - 1) * gridSize + surfCharge[0]*Psi[0] + surfCharge[1]*Psi[M-1];

    for (int j = 0; j < M; j++)
    {
        temp[j] = rhoZ[j]*Y[j];
    }

    Int3 = 0.5*trapz(temp, 0, M - 1)*gridSize;

    if (rhoB[0] < 0.499999)
    {
        mu_b = log(rhoB[0] / (1 - 2 * rhoB[0]));
    }
    else
    {
        mu_b = 0;
    }

    for (int j = 0; j < M; j++)
    {
        temp[j] = mu_b*(rho[0][j]+rho[1][j]);
    }

    Int4 = trapz(temp, 0, M - 1)*gridSize;

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < M; j++)
        {
            temp[j] = H[i][j] * rho[i][j];
        }

        Int5 += trapz(temp, 0, M - 1) * gridSize;
    }

    W = (Int1 + Int2 + Int3 - Int4 - Int5)/separation;
    W_entr = Int1/separation;
    W_elec = Int2/separation;
    W_yuka = Int3/separation;
    W_chem = Int4/separation;
    W_ads = Int5/separation;

    delete[] temp;
}

void ILFT::calcGrandPotentialUniform()
{
    double Int1, Int2, Int3, Int4, Int5=0;
    double *temp;
    double mu_b;

    temp = new double[M];

    // Entropy Integrand
    for (int j = 0; j < M; j++)
    {
        if (rhoB[0] < 0.499999)
        {
            temp[j] = rhoB[0] * log(rhoB[0]) + rhoB[1] * log(rhoB[1]) + (1 - rhoB[0] - rhoB[1]) * log(1 - rhoB[0] - rhoB[1]);
        }
        else
        {
            temp[j] = rhoB[0] * log(rhoB[0]) + rhoB[1] * log(rhoB[1]);
        }
    }

    // Calculate entropy contribution
    Int1 = trapz(temp, 0, M - 1) * gridSize;

    // Electrostatic Integrand
    for (int j = 0; j < M; j++)
    {
        temp[j] = 0;
    }

    Int2 = trapz(temp, 0, M - 1) * gridSize;

    for (int j = 0; j < M; j++)
    {
        temp[j] = 0;
    }

    Int3 = 0.5*trapz(temp, 0, M - 1) * gridSize;

    if (rhoB[0] < 0.4999999)
    {
        mu_b = log(rhoB[0] / (1 - 2 * rhoB[0]));
    }
    else
    {
        mu_b = 0;
    }

    for (int j = 0; j < M; j++)
    {
        temp[j] = mu_b * (rhoB[0] + rhoB[1]);
    }

    Int4 = trapz(temp, 0, M - 1) * gridSize;

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < M; j++)
        {
            temp[j] = H[i][j] * rho[i][j];
        }

        Int5 += trapz(temp, 0, M - 1) * gridSize;
    }

    Wuniform = (Int1 + Int2 + Int3 - Int4 - Int5)/separation;

    delete[] temp;
}