#ifndef ILFT_H
#define ILFT_H

/* Constants */
#define epsilon0 8.854187817E-12L
#define kB 1.3806504E-23L
#define e0 1.60217646E-19L
#define Na 6.02214179E+23L
#define Pi 3.1415926536

#include "../global_external.h"

/* Include utilities */
#include "../utilities/utilities.h"

class ILFT
{
private:
    /*=============================
    DFT Variables
    =============================*/
    /* System variables */
    int M;             // Total number of grid points
    int M_old;         // Total number of grid points on last iteration
    double lunit;      // Length unit
    double separation; // Starting separation (length units)
    double gridSize;   // Grid sizing (length units)
    double T;          // Temperature (K)
    double eps;        // Dielectric constant (relative)
    double BJ;         // Bjerrum length (length units)
    double lambdaD;	   // \lambda_D from the theory
    double kappa;      // Inverse Debye length
    double alpha;      // first constant for Yukawa
    double ell;        // correlation length for Yukawa (length units)
    double eta;        // Strength of interaction between total density and pref adsoprtion potential
    double chi;        // second constant for Yukawa
    double a;          // preferential adsorption range, usually the molecular size (b)
    double netcharge;  // Variable to hold the net charge for checking results

    /* Iterative variables */
    int maxIters;
    int numSteps;
    int itOuter;
    int iter;
    bool odd;
    bool init;
    double size_out;
    double dsize0;
    double dsize1;
    double errorTol;
    double error;
    double error_prev;
    bool BCswitch;
    double iterlow; // low value for outer iteration (e.g. phi, sigma, V iteration schemes)
    double iterhigh; // high value for outer iteration
    double numiter; // number of steps to do for outer iteration
    double iterstep; // amount to increment scanning variable, calculated from iterlow, iterhigh, and numiter
    short itertype;  // Sets what variable should be iterated over, if any
    double iterval; // current value of the parameter being scanned
    bool booldir; // Sets the directions that you are scanning when scanning phi (false=forward;True=backward)
    std::string dir; // name of directory that data will be output to when going either forward or reverse
    int count; // set counter so that when iterating over phi, you stop after forward+reverse and dont go forever 


    /* Energy Calculation */
    double Ener0;
    double Ener1;

    /* Poisson equation */
    double xPsi;

    double W;
    double W_entr;
    double W_elec;
    double W_yuka;
    double W_chem;
    double W_ads;
    double mu;
    double minW;
    double Wuniform;

    /* Mixing density */
    long double mixF;
    long double minMix;
    long double maxMix;
    long double maxMixBase;
    long double mixCoeff;

    /* Boundary variables */
    short BC_type;  // Type of boundary condition
    double Qmin;	// Minimum Q
    double Qmax;	// Maximum Q
    double Qn;		// Number of Q to do
    double Qstep;
    double *surfBC; // Surface boundary condition
    double *surfCharge;
    double *sigma_initial;
    double *sigma_target;
    
    /* Arrays for each component */
    short *Z;       // Valency
    double *h;
    double *D;      // Diameter (length unit)
    //double *mu;     // Total chemcial potential
    double *rhoB;   // Bulk density (reduced density)
    double **uSurf; // Interaction parameter with each surface

    /* Spatially defined arrays */
    double *Psi;
    double *PsiNew;
    double *Y;
    double *YNew;
    double *rhoZ;
    double *rhoZf;
    double **H;
    double **Ext;
    double **rho;
    double **rhof;
    double **rhoRef;
    double **rhoNew;
    double **rhoTemp;

    /* Conversion constants */
    double mol_to_red;
    double V_to_red;
    double nm2_to_red;

    /*=============================
      Data I/O
    =============================*/
    /* Directories */
    std::string cwd;       // current working directory
    std::string outputDir; // base level directory for data output (relative to cwd)

    /* File names */
    std::string logFilename = "log.out";
    std::string iterFilename = "iterations.out";
    std::string pressureDAT = "pressure.dat";
    std::string energyDAT = "energy.dat";
    std::string failedFilename = "failed.dat";
    std::string iterationDAT;
    std::string densityDAT;
    std::string densityInitDAT;
    std::string surfacechargeDAT;
    std::string freeEnergyDAT;
    std::string systemName;
    std::string systemNameFileName;

    /* Actual Files */
    std::string inputFile;
    std::string logFile;
    std::string iterFile;
    std::string pressureFile;
    std::string energyFile;
    std::string densityFile;
    std::string surfacechargeFile;
    std::string freeEnergyFile;
    std::string initFile;
    std::string failedFile;

    /* std::ofstream to do the actual I/O */
    std::ifstream inputStream;
    std::ifstream previousSystemNameStream;
    std::ifstream previousSystemStream;
    std::ofstream logStream;
    std::ofstream failedStream;
    std::ofstream iterStream;
    std::ofstream pressureStream;
    std::ofstream energyStream;
    std::ofstream densityStream;
    std::ofstream surfacechargeStream;
    std::ofstream freeEnergyStream;
    std::ofstream systemNameStream;

    void setOutputFileStructure();
    void closeFiles();
    void readInput(std::string inputPath);

    /*=============================
      Numerical tools
    =============================*/
    double trapz(double *f, int lower_bound, int upper_bound);
    double simps1(double *f, int lower_bound, int upper_bound);

public:
    /* Default constructor/destructor*/
    ILFT();
    ~ILFT();

    /* Constructor */
    ILFT(std::string inputPath);

    /*=============================
      Main Engine
    =============================*/
    void initializeSystem();
    void solveDensity();
    void finalizeSystem();

    void updateDensity();
    void mixDensity();
    void mixPotential();
    void mixY();
    void calculateError();
    void checkChargeNeutrality();

    /*=============================
      Calculations
    =============================*/
    void calcChargeDensity();
    void calcChargeDensityf();
    void weightedDensity(int i, double **NI);
    void calcChemicalPotential();
    void calcDensity(double x);
    void calcDensityInitial();
    void solveYukawaf();
    void calcSurfaceCharge();
    void calcGrandPotential();
    void calcGrandPotentialUniform();
    void calcAdsorptionProfiles();
    int getsgn(double x);
    void calcNetCharge();


    /*=============================
      Poisson Equation
    =============================*/
    void solvePoisson();
    void solvePoissonPotential();
    void solvePoissonSurfaceCharge();
    void solvePoissonSurfaceCharge2();
    void chargeNeutral();

    /*=============================
      Yukawa Field Equation
    =============================*/
    void solveYukawa();
    void solveYukawa2();

    /*=============================
      Iterative controls
    =============================*/
    void increaseSeparation();
    void increaseSurfaceCharge();
    void stepRhoB();
    void stepSC();
    void stepV();
    void getIterval();
    void transition();
    void flipdir();

    /*=============================
      Send data to main routine
    =============================*/
    int getSteps();

    /*=============================
      Data I/O
    =============================*/
    void writeInitial();
    void openLog();
    void writeLog(std::string logText);
    void openIter();
    void writeIter(int iter);
    void closeIter();
    void writeDensPot();
    void writeSurfaceCharge();
    void writeFreeEnergy();
    void writeDensPotInitial();
    void writeRhof();
    void finalizeIO();
    void readInitSystemName(std::string inputPath);
    void readInitSystemData();
    void writeProfileFilename();
};

#endif //ILFT_H