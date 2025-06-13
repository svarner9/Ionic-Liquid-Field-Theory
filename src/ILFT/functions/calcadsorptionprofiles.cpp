#include "../ILFT.h"

void ILFT::calcAdsorptionProfiles()
{

    int b = D[0];

    int left = ceil(b/gridSize);
    int right = M-(left+1);

    // std::cout << "Entered calcAdsorptionProfiles" << std::endl;


    for (int i=0; i<2; i++)
    {

        for (int j = 0; j <= left; j++)
        {
            H[i][j] = h[i] * pow(1.0 - (double)j / left, 2);
            H[i][M - (j + 1)] = h[i] * pow(1.0 - (double)j / left, 2);
        }

        // Profile in center
        for (int j = left + 1; j < right; j++)
        {
            H[i][j] = 0;
        }
    }

    // std::cout << "Exited calcAdsorptionProfiles" << std::endl;
}