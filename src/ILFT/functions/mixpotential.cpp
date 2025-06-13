#include "../ILFT.h"

void ILFT::mixPotential()
{
	/*=============================
       Calculate mixing parameter based on error
     =============================*/
	// if (error > error_prev)
	// {
	// 	mixF = minMix;
	// 	maxMix = maxMixBase * 0.1;
	// }
	// else if (error > 1.0E0)
	// {
	// 	mixF = minMix;
	// 	maxMix = maxMixBase * 0.1;
	// }
	// else if (error > 1.0E-1)
	// {
	// 	mixF += 0.1 * mixCoeff;
	// 	maxMix = maxMixBase;
	// }
	// else if (error > 1.0E-2)
	// {
	// 	mixF += 0.3 * mixCoeff;
	// 	maxMix = maxMixBase * 5;
	// }
	// else if (error > 1.0E-3)
	// {
	// 	mixF += 0.5 * mixCoeff;
	// 	maxMix = maxMixBase * 7;
	// }
	// else if (error > 1.0E-4)
	// {
	// 	mixF += 0.7 * mixCoeff;
	// 	maxMix = maxMixBase * 50;
	// }
	// else if (error > 1.0E-5)
	// {
	// 	mixF += 0.9 * mixCoeff;
	// 	maxMix = maxMixBase * 100;
	// }
	// else
	// {
	// 	mixF += 1 * mixCoeff;
	// 	maxMix = maxMixBase * 1000;
	// }

	// if (mixF > maxMix)
	// {
	// 	mixF = maxMix;
	// }

	// if (iter < 5E5)
	// {
	// 	mixF = minMix;
	// }

	/*=============================
      Picard Mixing
    =============================*/

	double asdf = 1E-7;

	for (int j = 0; j < M; j++)
	{
		Psi[j] = (1.0 - asdf) * Psi[j] + asdf * PsiNew[j];
	}

}

void ILFT::mixY()
{
	/*=============================
       Calculate mixing parameter based on error
     =============================*/
	// if (error > error_prev)
	// {
	// 	mixF = minMix;
	// 	maxMix = maxMixBase * 0.1;
	// }
	// else if (error > 1.0E0)
	// {
	// 	mixF = minMix;
	// 	maxMix = maxMixBase * 0.1;
	// }
	// else if (error > 1.0E-1)
	// {
	// 	mixF += 0.1 * mixCoeff;
	// 	maxMix = maxMixBase;
	// }
	// else if (error > 1.0E-2)
	// {
	// 	mixF += 0.3 * mixCoeff;
	// 	maxMix = maxMixBase * 5;
	// }
	// else if (error > 1.0E-3)
	// {
	// 	mixF += 0.5 * mixCoeff;
	// 	maxMix = maxMixBase * 7;
	// }
	// else if (error > 1.0E-4)
	// {
	// 	mixF += 0.7 * mixCoeff;
	// 	maxMix = maxMixBase * 50;
	// }
	// else if (error > 1.0E-5)
	// {
	// 	mixF += 0.9 * mixCoeff;
	// 	maxMix = maxMixBase * 100;
	// }
	// else
	// {
	// 	mixF += 1 * mixCoeff;
	// 	maxMix = maxMixBase * 1000;
	// }

	// if (mixF > maxMix)
	// {
	// 	mixF = maxMix;
	// }

	// if (iter < 5E5)
	// {
	// 	mixF = minMix;
	// }

	/*=============================
      Picard Mixing
    =============================*/

	double asdf = 1E-7;

	for (int j = 0; j < M; j++)
	{
		Y[j] = (1.0 - asdf) * Y[j] + asdf * YNew[j];
	}
}