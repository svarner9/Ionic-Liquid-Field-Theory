//********************************************************************//
// Ionic Liquid Field Theory                                          //
// Authors: Sam Varner / Chris Balzer (svarner@caltech.edu)           //
// Date: April 2021                                                   //
//********************************************************************//

#include "global_external.h"
#include "global_internal.h"

using namespace SysInfo;
using namespace std;

int main(int argc, char *argv[])
{
    // Get the number of OpenMP threads
    GET_OMP_THREADS();

    // Initialize timers
    Timers MainTimer;
    // Timers IterationTimer;

    // Parse input arguments
    string inputPath;
    if (argc > 1)
    {
        inputPath = argv[1];
    }
    else
    {
        inputPath = "./";
    }

    // Initialize DFT engine
    std::cout << "Initializing system..." << std::endl;
    ILFT ilft = ILFT(inputPath);

    // Start timer for main routine
    START_CPU_TIMER(MainTimer);

	
    

	// Solve density profile using engine
	ilft.solveDensity();

	// Stop iteration timer
	// STOP_CPU_TIMER(IterationTimer);

    // Average iteration time
    // IterationTimer.accum /= ilft.getSteps() + 1;

    // Stop timer for main routine
    STOP_CPU_TIMER(MainTimer);

    // Write timing to output
    int number1 = MainTimer.duration;
    std::stringstream ss1;
    ss1 << number1;
    ilft.writeLog("\nTotal Program Time   : " + ss1.str() + " seconds");
    
    int number2 = MainTimer.duration / (ilft.getSteps()+1);
    std::stringstream ss2;
    ss2 << number2;
    ilft.writeLog("Average Time per Run : " + ss2.str() + " seconds");

    // Finalize I/O
    ilft.finalizeIO();

    // Print the total time to the terminal as well
    std::cout << "\nTotal Program Time : " << MainTimer.duration << " seconds" << std::endl;

    return 0;
}