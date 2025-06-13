#include "../ILFT.h"

void ILFT::readInput(std::string inputPath)
{

    // Initialize some variables
    int item;
    std::string s;
    std::string delim = ";";
    std::string sub_s;
    std::size_t pos = 0;

    // Set the input file name
    inputFile = inputPath + "/input.dat";

    // Open the file
    inputStream.open(inputFile.c_str(), std::fstream::in);

    // Check if file has been opened
    if (inputStream)
    {
        // Read through input file
        while (std::getline(inputStream, s))
        {

            if (s.compare("SIZE:") == 0)
            {
                if (std::getline(inputStream, s))
                {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos)
                    {
                        sub_s = s.substr(0, pos);

                        if (item == 0)
                            separation = atof(sub_s.c_str());
                        if (item == 1)
                            gridSize = atof(sub_s.c_str());
                        s.erase(0, pos + delim.length());
                        item++;
                    }
                }
                else
                {
                    std::cerr << "DFT::read_input - Could not read in SIZE" << std::endl;
                    std::exit(0);
                }
            }

            if (s.compare("BULK_DENSITY:") == 0)
            {
                if (std::getline(inputStream, s))
                {
                    // Store the value
                    rhoB[0] = atof(s.c_str());
                    rhoB[1] = atof(s.c_str());
                }
                else
                {
                    std::cerr << "DFT::read_input - Could not read in BULK_DENSITY" << std::endl;
                    std::exit(0);
                }
            }

            if (s.compare("SURFACE:") == 0)
            {
                if (std::getline(inputStream, s))
                {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos)
                    {
                        sub_s = s.substr(0, pos);

                        if (item == 0)
                            surfBC[0] = atof(sub_s.c_str());
                        if (item == 1)
                            surfBC[1] = atof(sub_s.c_str());
                        if (item == 2)
                            BC_type = atof(sub_s.c_str());
                        if (item == 3)
                            uSurf[0][0] = atof(sub_s.c_str());
                        if (item == 4)
                            uSurf[0][1] = atof(sub_s.c_str());
                        if (item == 5)
                            uSurf[1][0] = atof(sub_s.c_str());
                        if (item == 6)
                            uSurf[1][1] = atof(sub_s.c_str());

                        s.erase(0, pos + delim.length());
                        item++;
                    }
                }
                else
                {
                    std::cerr << "DFT::read_input - Could not read in SURFACE" << std::endl;
                    std::exit(0);
                }
            }

            if (s.compare("VALENCY:") == 0)
            {
                if (std::getline(inputStream, s))
                {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos)
                    {
                        sub_s = s.substr(0, pos);

                        if (item == 0)
                            Z[0] = atof(sub_s.c_str());
                        if (item == 1)
                            Z[1] = atof(sub_s.c_str());

                        s.erase(0, pos + delim.length());
                        item++;
                    }
                }
                else
                {
                    std::cerr << "DFT::read_input - Could not read in VALENCY" << std::endl;
                    std::exit(0);
                }
            }

            if (s.compare("DIAMETER:") == 0)
            {
                if (std::getline(inputStream, s))
                {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos)
                    {
                        sub_s = s.substr(0, pos);

                        if (item == 0)
                            D[0] = atof(sub_s.c_str());
                        if (item == 1)
                            D[1] = atof(sub_s.c_str());

                        s.erase(0, pos + delim.length());
                        item++;
                    }
                }
                else
                {
                    std::cerr << "DFT::read_input - Could not read in DIAMETER" << std::endl;
                    std::exit(0);
                }
            }

            if (s.compare("SYSTEM:") == 0)
            {
                if (std::getline(inputStream, s))
                {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos)
                    {
                        sub_s = s.substr(0, pos);

                        if (item == 0)
                            T = atof(sub_s.c_str());
                        if (item == 1)
                            lunit = atof(sub_s.c_str());
                        if (item == 2)
                            eps = atof(sub_s.c_str());
                        if (item == 3)
                            alpha = atof(sub_s.c_str());
                        if (item == 4)
                            ell = atof(sub_s.c_str());
                        if (item == 5)
                            h[0] = atof(sub_s.c_str());
                        if (item == 6)
                            h[1] = atof(sub_s.c_str());

                        s.erase(0, pos + delim.length());
                        item++;
                    }
                }
                else
                {
                    std::cerr << "DFT::read_input - Could not read in SYSTEM" << std::endl;
                    std::exit(0);
                }
            }

            if (s.compare("ITERATIVE:") == 0)
            {
                if (std::getline(inputStream, s))
                {
                    item = 0;
                    pos = 0;
                    while ((pos = s.find(delim)) != std::string::npos)
                    {
                        sub_s = s.substr(0, pos);
                        if (item == 0)
                            maxIters = atof(sub_s.c_str());
                        if (item == 1)
                            errorTol = atof(sub_s.c_str());
                        if (item == 2)
                            itertype = atof(sub_s.c_str());
                        if (item == 3)
                            numiter = atof(sub_s.c_str());
                        if (item == 4)
                            iterlow = atof(sub_s.c_str());
                        if (item == 5)
                            iterhigh = atof(sub_s.c_str());

                        s.erase(0, pos + delim.length());
                        item++;
                    }
                }
                else
                {
                    std::cerr << "DFT::read_input - Could not read in ITERATIVE" << std::endl;
                    std::exit(0);
                }
            }

            if (s.compare("OUT_PATH:") == 0)
            {
                if (std::getline(inputStream, s))
                {
                    // Store the value
                    outputDir = outputDir + s;
					// std::cout << " ___ " << outputDir << " \n";
                }
                else
                {
                    std::cerr << "DFT::read_input - Could not read in OUT_PATH" << std::endl;
                    std::exit(0);
                }
            }
        }

        // Close the input file
        inputStream.close();
    }
    else
    {
        std::cout << "DFT::read_input - No input file found. Using default." << std::endl;
    }

    // Set up the file outputs
    // setOutputFileStructure();
    // openLog();
}

void ILFT::writeInitial()
{
    // Write initial system
    logStream << "=======================DFT for Electrolyte System========================"
              << " \n\n";
    logStream << "              Copyright Owners: Chris Balzer (May 2020)"
              << " \n";
    logStream << "========================================================================="
              << " \n";

    logStream << "\n==============="
              << " \n";
    logStream << "INITIALIZATION"
              << " \n";
    logStream << "==============="
              << " \n";
    logStream << "Data read from input file   : " << inputFile << " \n";
    logStream << "Writing output to directory : " << outputDir << " \n";

    /* System parameters from input file */
    logStream << "    Separation     		  : " << separation << " \n";
    logStream << "    Ion Bulk Vol Frac	      : " << rhoB[0] << " \n";
    logStream << "    Valence  (Pos : Neg)    : " << Z[0] << " : " << Z[1] << " \n";
    logStream << "    Diameter (Pos : Neg)    : " << D[0] << " : " << D[1] << " \n";
    logStream << "    Iteration type/var      : " << itertype << " \n";
    logStream << "    Starting iter value     : " << iterlow << " \n";
	logStream << "    Ending iter value       : " << iterhigh << " \n";
	logStream << "    Number of steps         : " << numiter << " \n";
    logStream << "    Length Unit             : " << lunit << " \n";
    logStream << "    Grid Sizing             : " << gridSize << " \n";
    logStream << "    Temperature (K)         : " << T << " \n";
    logStream << "    Dielectric Constant     : " << eps << " \n";
    logStream << "    Alpha                   : " << alpha << " \n";
    logStream << "    Ell_0                   : " << ell << " \n";
    logStream << "    h_a, h_c                : " << h[0] << "," << h[1] << " \n";
    logStream << "    Error Tolerance         : " << errorTol << " \n";

    logStream << "\n==============="
              << " \n";
    logStream << "RUNNING ILFT"
              << " \n";
    logStream << "==============="
              << std::endl;
}

void ILFT::writeLog(std::string logText)
{
    logStream << logText << " \n";
}

void ILFT::openIter()
{
	std::stringstream ss3;
	ss3 << iterval;

	iterFile = outputDir + dir + iterationDAT + ss3.str() + ".out";

	iterStream.open(iterFile.c_str(), std::fstream::out | std::ofstream::trunc);
}

void ILFT::writeIter(int iter)
{
    iterStream << "Iteration: " << iter << "  ";
    iterStream << "Error: " << std::scientific << std::setprecision(4) << error << "  ";
    iterStream << "Mix Coeff: " << std::scientific << std::setprecision(4) << mixF << "  ";
    iterStream << "Free Energy: " << std::scientific << std::setprecision(6) << W << std::endl;
}

void ILFT::closeIter()
{
    iterStream.close();
}

void ILFT::writeDensPot()
{
    // Get file name and open file
    std::stringstream ss3;
    ss3 << iterval;

    densityFile = outputDir + dir + densityDAT + ss3.str() + ".dat";
    densityStream.open(densityFile.c_str(), std::fstream::out | std::ofstream::trunc);

    // Print to file
    for (int j = 0; j < M; j++)
    {
        densityStream << std::fixed << std::setprecision(5) << j * gridSize * lunit / 1.0E-9 << " ";
        densityStream << std::fixed << std::setprecision(8) << rho[0][j] << " ";
        densityStream << std::fixed << std::setprecision(8) << rho[1][j] << " ";
        densityStream << std::fixed << std::setprecision(8) << Psi[j] << " ";
        densityStream << std::fixed << std::setprecision(8) << Y[j] << " ";
        densityStream << std::fixed << std::setprecision(5) << H[0][j] << " ";
        densityStream << std::fixed << std::setprecision(5) << H[1][j] << std::endl;
    }

    // Close the file
    densityStream.close();
}

void ILFT::writeRhof()
{
    failedFile = outputDir + failedFilename;
    failedStream.open(failedFile.c_str(), std::fstream::out | std::ofstream::trunc);

    for (int j=0; j < M; j++)
    {
        failedStream << std:: fixed << std::setprecision(5) << j * gridSize * lunit / 1.0E-9 << " ";
        failedStream << std:: fixed << std::setprecision(8) << rhof[0][j] << " ";
        failedStream << std:: fixed << std::setprecision(8) << rhof[1][j] << " ";
        failedStream << std::fixed << std::setprecision(8) << Psi[j] << " ";
        failedStream << std::fixed << std::setprecision(8) << Y[j] << std::endl;

    }

    failedStream.close();
}

void ILFT::writeSurfaceCharge()
{
    // Get file name and open file
    std::stringstream ss3;
    ss3 << iterval;

    surfacechargeFile = outputDir + dir + surfacechargeDAT + ss3.str() + ".dat";
    surfacechargeStream.open(surfacechargeFile.c_str(), std::fstream::out | std::ofstream::trunc);

    // Print to file
    surfacechargeStream << std::fixed << std::setprecision(12) << surfCharge[0] << " ";
    surfacechargeStream << std::fixed << std::setprecision(12) << surfCharge[1] << std::endl;

    // Close the file
    surfacechargeStream.close();

}

void ILFT::writeFreeEnergy()
{
    // Get file name and open file
    std::stringstream ss3;
    ss3 << iterval;

    freeEnergyFile = outputDir + dir + freeEnergyDAT + ss3.str() + ".dat";
    freeEnergyStream.open(freeEnergyFile.c_str(), std::fstream::out | std::ofstream::trunc);

    // Print to file
    freeEnergyStream << std::scientific << std::setprecision(12) << mu << " ";
    // freeEnergyStream << std::scientific << std::setprecision(8) << W_entr << " ";
    // freeEnergyStream << std::scientific << std::setprecision(8) << W_elec << " ";
    // freeEnergyStream << std::scientific << std::setprecision(8) << W_yuka << " ";
    // freeEnergyStream << std::scientific << std::setprecision(8) << W_chem << " ";
    freeEnergyStream << std::scientific << std::setprecision(12) << W << " ";
    // freeEnergyStream << std::scientific << std::setprecision(8) << Wuniform << " \n";

    // Close the file
    freeEnergyStream.close();

}

void ILFT::writeDensPotInitial()
{
    // Get file name and open file
    std::stringstream ss3;
	std::stringstream ss4;
    ss3 << "init";
	ss4 << iterval;

    densityFile = outputDir + dir + densityInitDAT + ss4.str() + "_" + ss3.str() + ".dat";
    densityStream.open(densityFile.c_str(), std::fstream::out | std::ofstream::trunc);

    // Print to file
    for (int j = 0; j < M; j++)
    {
        densityStream << std::fixed << std::setprecision(5) << j * gridSize * lunit / 1.0E-9 << " ";
        densityStream << std::fixed << std::setprecision(8) << rho[0][j] << " ";
        densityStream << std::fixed << std::setprecision(8) << rho[1][j] << " ";
        densityStream << std::fixed << std::setprecision(8) << Psi[j] << " ";
        densityStream << std::fixed << std::setprecision(8) << Y[j] << " ";
        densityStream << std::fixed << std::setprecision(5) << H[0][j] << " ";
        densityStream << std::fixed << std::setprecision(5) << H[1][j] << std::endl;
    }

    // Close the file
    densityStream.close();
}

void ILFT::finalizeIO()
{
    // Write a message to end
    logStream << "\n==============="
              << " \n";
    logStream << "END"
              << " \n";
    logStream << "==============="
              << std::endl;

    // Close all appending files
    closeFiles();
}

void ILFT::closeFiles()
{
    logStream.close();
}

void ILFT::setOutputFileStructure()
{
    // Set file streams
	iterationDAT = "iterations/iterations_";
    densityDAT = "density_potential/density_potential_";
	densityInitDAT = "initial/density_potential_";
    surfacechargeDAT = "surface_charge/surface_charge_";
    freeEnergyDAT = "free_energy/free_energy_";

    // Make the output directory if it isn't already made
    struct stat sb;
    std::string outputBase = outputDir;
    std::string outputDirection = outputDir + dir;
    std::string outputDensPot = outputDir + dir + "density_potential/";
    std::string outputSurfCha = outputDir + dir + "surface_charge/";
    std::string outputInitial = outputDir + dir + "initial/";
    std::string outputIterations = outputDir + dir + "iterations/";
    std::string outputFreeEnergy = outputDir + dir + "free_energy/";

    if (stat(outputBase.c_str(), &sb) == -1)
    {
        std::string mkdirBase = "mkdir " + outputBase;
        int systemRet = std::system(mkdirBase.c_str());
        if (systemRet == -1)
        {
            std::cerr << "DFT : failed to create output directory" << std::endl;
            exit(0);
        }
    }

    if (stat(outputDir.c_str(), &sb) == -1)
    {
        std::string mkdirBase = "mkdir " + outputDir;
        int systemRet = std::system(mkdirBase.c_str());
        if (systemRet == -1)
        {
            std::cerr << "DFT : failed to create output directory" << std::endl;
            exit(0);
        }
    }

    if (stat(outputDirection.c_str(), &sb) == -1)
    {
        std::string mkdirBase = "mkdir " + outputDirection;
        int systemRet = std::system(mkdirBase.c_str());
        if (systemRet == -1)
        {
            std::cerr << "DFT : failed to create output + direction directory" << std::endl;
            exit(0);
        }
    }

    if (stat(outputDensPot.c_str(), &sb) == -1)
	{
		std::string mkdirBase = "mkdir " + outputDensPot;
		int systemRet = std::system(mkdirBase.c_str());
		if (systemRet == -1)
		{
			std::cerr << "DFT : failed to create density_potential directory" << std::endl;
			exit(0);
		}
	}

	if (stat(outputSurfCha.c_str(), &sb) == -1)
	{
		std::string mkdirBase = "mkdir " + outputSurfCha;
		int systemRet = std::system(mkdirBase.c_str());
		if (systemRet == -1)
		{
			std::cerr << "DFT : failed to create surface_charge directory" << std::endl;
			exit(0);
		}
	}

	if (stat(outputInitial.c_str(), &sb) == -1)
	{
		std::string mkdirBase = "mkdir " + outputInitial;
		int systemRet = std::system(mkdirBase.c_str());
		if (systemRet == -1)
		{
			std::cerr << "DFT : failed to create initial directory" << std::endl;
			exit(0);
		}
	}

	if (stat(outputIterations.c_str(), &sb) == -1)
	{
		std::string mkdirBase = "mkdir " + outputIterations;
		int systemRet = std::system(mkdirBase.c_str());
		if (systemRet == -1)
		{
			std::cerr << "DFT : failed to create iterations directory" << std::endl;
			exit(0);
		}
	}

	if (stat(outputFreeEnergy.c_str(), &sb) == -1)
	{
		std::string mkdirBase = "mkdir " + outputFreeEnergy;
		int systemRet = std::system(mkdirBase.c_str());
		if (systemRet == -1)
		{
			std::cerr << "DFT : failed to create free_energy directory" << std::endl;
			exit(0);
		}
	}

	
}

void ILFT::openLog()
{
    // Open log file
    logFile = outputDir + logFilename;
    logStream.open(logFile.c_str(), std::fstream::out | std::ofstream::trunc);
}

void ILFT::readInitSystemName(std::string inputPath)
{
    struct stat sb;
    std::string name_holder = inputPath + "/init_src.dat";
    std::string s;
    std::string filepath;

    std::cout << "name_holder = " << name_holder.c_str() << std::endl;

    if (stat(name_holder.c_str(), &sb) != -1)
    {
        init = true;
        std::cout << "init_src file found!" << std::endl;

        // Open the file
        previousSystemNameStream.open(name_holder.c_str(), std::fstream::in);

        // Check if file has been opened
        if (previousSystemNameStream)
        {
            // Read through input file
            std::getline(previousSystemNameStream, systemName);
            std::cout << systemName << std::endl;
        }
        else
        {
            init = false;
            std::cout << "init_src file not opened" << std::endl;
        }
    }
    else
    {
        std::cout << "No init_src file!" << std::endl;
    }
    
}

void ILFT::readInitSystemData()
{
    // std::cout << "Reading output from ... " << systemName << std::endl;

    initFile = "/extra/svarner/RTILSC/Yukawa_ads_nofield/output/" + systemName + "/density_potential_0.dat";
    // initFile = "../init_source_profiles.dat";
    std::string delim = " ";
    std::string s;
    int item;
    int i = 0;
    std::size_t pos = 0;
    std::string sub_s;

    std::cout << initFile << std::endl;

    //Open the output file of the previous run
    previousSystemStream.open(initFile.c_str(), std::fstream::in);

    if (previousSystemStream)
    {
        std::cout << "Found and opened init file..." << std::endl;

        while (std::getline(previousSystemStream, s))
        {
            item = 0;
            pos = 0;
            while ((pos = s.find(delim)) != std::string::npos)
            {
                sub_s = s.substr(0, pos);

                if (item == 1)
                {
                    rho[0][i] = atof(sub_s.c_str());
                }
                if (item == 2)
                {
                    rho[1][i] = atof(sub_s.c_str());
                }
                if (item == 3)
                {
                    Psi[i] = atof(sub_s.c_str());
                }
                if (item == 4)
                {
                    Y[i] = atof(sub_s.c_str());
                }
                if (item == 5)
                {
                    H[0][i] = atof(sub_s.c_str());
                }
                if (item ==  6)
                {
                    H[1][i] = atof(sub_s.c_str());
                }
                s.erase(0, pos + delim.length());
                item++;
            }
            i++;
        }
    }
    else
    {
        std::cout << "Not able to open init file..." << std::endl;
    }
}