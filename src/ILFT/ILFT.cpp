#include "ILFT.h"

/* Default constructor */
ILFT::ILFT() = default;

/* Default destructor */
ILFT::~ILFT()
{

    for (int i = 0; i < 2; i++)
    {
        delete[] uSurf[i];
    }

    delete[] Z;
    delete[] rhoB;
    delete[] D;
    delete[] surfBC;
    delete[] uSurf;
    delete[] h;
}

/* Constructor */
ILFT::ILFT(std::string inputPath)
{
    // Initialize variables
    Z = new short[2];
    rhoB = new double[2];
    h = new double[2];
    D = new double[2];
    surfBC = new double[2];

    uSurf = new double *[2];
    for (int i = 0; i < 2; i++)
    {
        uSurf[i] = new double[2];
    }

    outputDir = "";

    //  Read input file and open files
    readInput(inputPath);

    if (itertype == 0)
    {
        count = 1;
        if (rhoB[0] > 0.49)
        {
            booldir = true; // sets direction to reverse
            dir = "reverse/";
        }
        else
        {
            booldir = false; // sets direction to forward
            dir = "forward/";
        }
    }
    else
    {
        dir = "";
        count = 2;
    }

    // Calculate gridpoints and sizing
    M = ceil(separation / gridSize) + 1;
    gridSize = separation / ((double)M - 1);

    // Set the iterval which is used to create output filenames
    getIterval();

    setOutputFileStructure();
    openLog();

    // Write initial message with input parameters
    writeInitial();

    // Calculate conversion factors and convert inputs
	// I changed these to 1 so I could just make inputs in reduced units
    mol_to_red = 1; // 1.0E3 * lunit * lunit * lunit * Na;
    V_to_red = 1;//e0 / (kB * T);
    nm2_to_red = lunit * lunit / (1.0E-9 * 1.0E-9);

    for (int i = 0; i < 2; i++)
    {
        rhoB[i] *= mol_to_red;
    }

    calcChemicalPotential();

    // if (BC_type == 0)
    // {
    //     for (int i = 0; i < 2; i++)
    //     {
    //         surfBC[i] *= V_to_red;
    //     }
    // }
    // else
    // {
    //     for (int i = 0; i < 2; i++)
    //     {
    //         surfBC[i] *= nm2_to_red;
    //     }
    // }

    std::cout << "surfBC[0] = " << surfBC[0] << std::endl;
    std::cout << "surfBC[1] = " << surfBC[1] << std::endl;

    // Calculate Bjerrum length and Debye length
    BJ = 0.25 * e0 * e0 / (Pi * lunit * eps * epsilon0 * kB * T);
	lambdaD = pow((4 * Pi * BJ),-0.5);
    kappa = sqrt(4 * Pi * BJ * (Z[0] * Z[0] * rhoB[0] + Z[1] * Z[1] * rhoB[1]));

    // Other variables
	iterstep = (iterhigh - iterlow) / numiter;
	odd = false;
    BCswitch = false;
    itOuter = 0;
	iter = 0;
    size_out = separation * lunit / 1.0E-9;
    mixF = 1E-4;
    minMix = 1E-4;
    maxMix = 1E-3;
    maxMixBase = 1e-3;
    mixCoeff = 1E-4;
    dsize0 = 2.0 * gridSize;
    dsize1 = 20.0 * gridSize;
    error = 0.01;
    xPsi = 0.0;
}

/* Main engine for DFT */
void ILFT::solveDensity()
{
    // Initialize run
    initializeSystem();

	while (itOuter <= numiter && count < 3)
	{

		// Initialize variables
		SysInfo::Timers timer;

		// Start timer
		START_CPU_TIMER(timer);

		// Open the iteration file for this surface charge
		openIter();

		// Write initial state to an output file
		writeDensPotInitial();

		// Calculate initial surface charge
		calcSurfaceCharge();
		
		// Update density until convergence
		do
		{
			// Calculate new density
			updateDensity();

			// Calculate error
			calculateError();

			// Mix densities
			mixDensity();

            // Calculate charge density
            calcChargeDensity();

			// Calculate the free energy of the system
			calcGrandPotential();

			// Output iteration info to iterative file
			if (iter % 10000 == 0)
			{
				writeIter(iter);
			}

			if (iter % 100000 == 0)
			{
				//Output intermediate density profile to file
				writeDensPot();

				calcSurfaceCharge();

				//Output surface charge to file
				writeSurfaceCharge();
			}

			// Increase iteration number
			iter++;
		} while (error > errorTol && iter < maxIters);

		// Stop timer
		STOP_CPU_TIMER(timer);

        calcChargeDensity();

        calcNetCharge();

		// Write to log file
        logStream << "Complete : Iterval = " << std::fixed << std::setprecision(5) << iterval << ", net charge = " << std::setprecision(5) << netcharge << ", duration (s) = " << std::setprecision(3) << timer.duration << " \n";

        transition();
    }

    // Finalize run
    finalizeSystem();
}

void ILFT::increaseSeparation()
{
    if (!odd)
    {
        separation += dsize1;
    }
    else
    {
        separation += dsize0;
        size_out = separation * lunit / 1.0E-9;
    }
}


void ILFT::transition()
{

	// Close the iteration file for current Q
	writeIter(iter);
    closeIter();

    // Calculate the Surface Charge
	calcSurfaceCharge();

	// Calculate the free energy for a uniform profile at the current bulk density
	calcGrandPotentialUniform();
    calcChemicalPotential();

	// Write the density and potentials
	writeDensPot();

	// Write the Surface Charge
	writeSurfaceCharge();

	// Write free energy to file
	writeFreeEnergy();

    if (itOuter == numiter && itertype == 0)
    {
        flipdir();
    }
    else
    {
        itOuter++;
        
        if (itertype == 0)
        {
            stepRhoB();
        }
        else if (itertype == 1)
        {
            stepSC();
        }
        else if (itertype == 2)
        {
            stepV();
        }
        else
        {
            std::cout << "Error with stepping in transition"
                      << " \n";
        }
    }

    getIterval();

    mixF = 1E-4;
    minMix = 1E-4;
    maxMix = 1E-3;
    maxMixBase = 1e-3;
    mixCoeff = 1E-4;
    error = 0.01;
	xPsi = 0.0;
	iter = 0;
    

    double t = 1.0e-6;
    // Perturb potentials a little bit to keep them nonzero
    for (int j = 0; j < M; j++)
    {
        Psi[j] += t - 2.0 * j / (M - 1) * t;
        PsiNew[j] += t - 2.0 * j / (M - 1) * t;
        Y[j] += t - 2.0 * j / (M - 1) * t;
        YNew[j] += t - 2.0 * j / (M - 1) * t;
    }

    calcDensityInitial();
    // Calculate the charge density
    calcChargeDensity();
}

void ILFT::initializeSystem()
{
    /* Calculate number of grid points */
    Psi = new double[M];
	PsiNew = new double[M];
    Y = new double[M];
	YNew = new double[M];
    rhoZ = new double[M];
    rhoZf = new double[M];
    Ext = new double *[2];
    rho = new double *[2];
    rhof = new double *[2];
    H = new double *[2];
    rhoRef = new double *[2];
    rhoNew = new double *[2];
    surfCharge = new double[2];

    minW = 1e10;
    W = 1e10;

    for (int i = 0; i < 2; i++)
    {
        Ext[i] = new double[M];
        rho[i] = new double[M];
        rhof[i] = new double [M];
        H[i] = new double[M];
        rhoRef[i] = new double[M];
        rhoNew[i] = new double[M];

        if (BC_type == 0)
        {
            for (int i = 0; i < 2; i++)
            {
                surfCharge[i] = 0.0;
            }
        }
        else
        {
            for (int i = 0; i < 2; i++)
            {
                surfCharge[i] = surfBC[i];
            }
        }
    }

    double t = 1.0e-6;

    // Initialize potentials
    for (int j = 0; j < M; j++)
    {
		Psi[j] = t - 2.0 * j  / (M - 1) * t;
        PsiNew[j] = t - 2.0 * j / (M - 1) * t;
        Y[j] = t - 2.0 * j / (M - 1) * t;
        YNew[j] = t - 2.0 * j / (M - 1) * t;
        rhoZ[j] = 0.0;
        rhoZf[j] = 0.0;
    }
	
	// Initialize density profiles to be uniform for now
    // for (int i = 0; i < 2; i++)
    // {
    for (int j = 0; j < M/2; j++)
    {
        rho[0][j] = rhoB[0];
        rho[1][j] = 0.0;
    }
    for (int j = M/2; j < M; j++)
    {
        rho[0][j] = 0.0;
        rho[1][j] = rhoB[1];
    }
    // }

    // calcDensityInitial();

	// Calculate the charge density
	calcChargeDensity();

    calcAdsorptionProfiles();

	// Calculate the external potential
	// Not currently running with any external potential
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < M; j++)
        {
            Ext[i][j] = 0.0;
        }
    }

}

void ILFT::finalizeSystem()
{
    /* Delete variables */
    for (int i = 0; i < 2; i++)
    {
        delete[] Ext[i];
        delete[] rho[i];
        delete[] rhof[i];
        delete[] rhoNew[i];
        delete[] H[i];
    }
    delete[] rhoZ;
    delete[] rhoZf;
    delete[] Psi;
	delete[] PsiNew;
    delete[] Y;
	delete[] YNew;
    delete[] Ext;
    delete[] rho;
    delete[] rhof;
    delete[] rhoNew;
    delete[] surfCharge;

}

int ILFT::getSteps()
{
    return numiter;
}