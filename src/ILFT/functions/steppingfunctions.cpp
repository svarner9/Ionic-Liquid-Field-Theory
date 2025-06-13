#include "../ILFT.h"

void ILFT::stepRhoB()
{
    if (booldir)
    {
        for (int i = 0; i < 2; i++)
        {
            rhoB[i] -= iterstep;
        }
    }
    else
    {
        for (int i = 0; i < 2; i++)
        {
            rhoB[i] += iterstep;
        }
    }

}

void ILFT::stepSC()
{
    for (int i = 0; i < 2; i++)
    {
        if (getsgn(surfBC[i]) == 0)
        {
            surfBC[i] += (i == 0) * iterstep - (i == 1) * iterstep;
        }
        else
        {
            surfBC[i] += getsgn(surfBC[i]) * iterstep;
        }
        
    }
}

void ILFT::stepV()
{
    for (int i = 0; i < 2; i++)
    {
        if (getsgn(surfBC[i]) == 0)
        {
            surfBC[i] += (i == 0) * iterstep - (i == 1) * iterstep;
        }
        else
        {
            surfBC[i] += getsgn(surfBC[i]) * iterstep;
        }
    }
}

void ILFT::getIterval()
{
    // Set what iterval is based on the itertype
    //  itertype = 0  --> iterval is rhoB
    //  itertype = 1 -->  iterval is SC
    //  itertype = 2  --> iterval is V
    if (itertype == 0)
    {
        iterval = rhoB[0];
    }
    else if (itertype == 1)
    {
        iterval = surfBC[0];
    }
    else if (itertype == 2)
    {
        iterval = surfBC[0];
    }
    else
    {
        std::cout << "Error with the iterval definition" << " \n";
    }
}

void ILFT::flipdir()
{

    if (booldir)
    {
        booldir = false;
        dir = "forward/";
    }
    else
    {
        booldir = true;
        dir = "reverse/";
    }

    itOuter = 0;
    count++;
    // Create new directories for other direction if not already created
    setOutputFileStructure();
}