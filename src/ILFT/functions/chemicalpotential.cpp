#include "../ILFT.h"

void ILFT::calcChemicalPotential()
{
	mu = log(rhoB[0] / (1 - 2 * rhoB[0]));
}