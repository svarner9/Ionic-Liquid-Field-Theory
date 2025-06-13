#include "../ILFT.h"

void ILFT::mixDensity()
{
    /*=============================
       Calculate mixing parameter based on error
     =============================*/
    if (error > error_prev)
    {
        mixF = minMix;
		maxMix = maxMixBase*0.5;
    }
	else if (error > 1.0E0)
	{
		mixF = minMix;
		maxMix = maxMixBase*0.5;
	}
    else if (error > 1.0E-1)
    {
        mixF += 0.1 * mixCoeff;
		maxMix = maxMixBase;
	}
    else if (error > 1.0E-2)
    {
        mixF += 0.3 * mixCoeff;
		maxMix = maxMixBase*5;
    }
    else if (error > 1.0E-3)
    {
        mixF += 0.5 * mixCoeff;
		maxMix = maxMixBase*7;
    }
    else if (error > 1.0E-4)
    {
        mixF += 0.7 * mixCoeff;
		maxMix = maxMixBase * 10;
    }
    else if (error > 1.0E-5)
    {
        mixF += 0.9 * mixCoeff;
		maxMix = maxMixBase * 10;
	}
	else if (error > 1.0E-6)
	{
		mixF += 1 * mixCoeff;
		maxMix = maxMixBase * 20;
	}
    else
    {
        mixF += 1 * mixCoeff;
		maxMix = maxMixBase * 25;
	}

    if (mixF > maxMix)
    {
        mixF = maxMix;
    }

	if (iter < 3E5)
	{
		mixF = minMix;
	}

    /*=============================
      Picard Mixing
    =============================*/
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < M; j++)
        {
			rho[i][j] = (1.0 - mixF) * rho[i][j] + mixF * rhoNew[i][j];
		}
    }

}