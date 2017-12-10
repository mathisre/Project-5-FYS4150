#include "math/random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using namespace std;
ofstream ofile;

int main(int numberOfArguments, char **argumentList)
{
    if (numberOfArguments != 3) {
        cout << "This program uses command line arguments to set the temperature and number of steps." << endl << "Set initial temperature in Kelvin and number of steps" << endl;
        exit(EXIT_FAILURE);
    }

    double initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[1]));
    int totalNumberOfSteps = atoi(argumentList[2]);

    double dt = UnitConverter::timeFromSI(1e-15);

    int numberOfUnitCells = 5;
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;

    System system;
    StatisticsSampler statisticsSampler;
    system.potential().setEpsilon(UnitConverter::energyFromSI(119.8*UnitConverter::kb));
    system.potential().setSigma(3.405);
    system.setPeriodicSize((numberOfUnitCells)*latticeConstant);

    system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
    system.removeTotalMomentum();

    double Initial_energy;
    double tol_E;
    double DiffConst;
    bool firststep;

    int timestep = 0;
    firststep = true;

    string filenameMovie = "argon_" + to_string(int(UnitConverter::temperatureToSI(initialTemperature))) + ".xyz";
    IO movie(filenameMovie);

    cout << endl << endl <<  setw(20) << "Timestep" <<
            setw(20) << "Time" <<
            setw(20) << "Temperature [K]" <<
            setw(20) << "KineticEnergy" <<
            setw(20) << "PotentialEnergy" <<
            setw(20) << "TotalEnergy" << endl;

    while (timestep < totalNumberOfSteps){
        system.step(dt);
        statisticsSampler.sample(system);

        if(firststep) {
            Initial_energy = statisticsSampler.totalEnergy();
            tol_E = 0.1;
            firststep = false;
        }

        if( timestep % 100 == 0 ) {
                cout << setw(20) << system.steps() <<
                        setw(20) << system.time() <<
                        setw(20) << statisticsSampler.temperature() <<
                        setw(20) << statisticsSampler.kineticEnergy() <<
                        setw(20) << statisticsSampler.potentialEnergy() <<
                        setw(20) << statisticsSampler.totalEnergy() << endl;
        }

        if( timestep % 5 == 0 ){ // Unit tests
            system.checkMomentum();
            if (statisticsSampler.totalEnergy() > tol_E + Initial_energy || statisticsSampler.totalEnergy() < Initial_energy - tol_E) {
                cout << "Error: energy not conserved after " << timestep << " timesteps." <<endl;
                cout << "Energy change  = " << statisticsSampler.totalEnergy() - Initial_energy << endl;
                exit(EXIT_FAILURE);
            }
        }
        movie.saveState(system);
        timestep++;
        } // End timestep loop
    movie.close();

    return 0;
}
