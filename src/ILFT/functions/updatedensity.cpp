#include "../ILFT.h"

void ILFT::updateDensity()
{

    /*=============================
       Calculate based on current density profile
     =============================*/
    // Solve for electric potential

	calcChargeDensity();

	solvePoisson();

	solveYukawa();

    /*=============================
      Substitute into mean field solutions for ion density
    =============================*/
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < M; j++)
        {
          rhoNew[i][j] = exp(-Z[i] * (Psi[j] + Y[j])) / ((1 / rhoB[i] - 2) * exp(-1 * H[i][j]) + 2 * cosh((Psi[j] + Y[j])));
        }
    }

    /*=============================
      Ensuring charge neutrality for fixed surface charge
    =============================*/
    // if (BC_type == 1)
    // {
    //   chargeNeutral();
    // }
}