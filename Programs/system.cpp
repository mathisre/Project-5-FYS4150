#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::RemoveLattice(){
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
    m_time = 0;
    m_steps = 0;
}

void System::applyPeriodicBoundaryConditions(){
    for(Atom *atom : m_atoms) {
        if (atom->position[0] > m_systemSize[0])  atom->position[0] -= m_systemSize[0];

        if (atom->position[1] > m_systemSize[1])  atom->position[1] -= m_systemSize[1];

        if (atom->position[2] > m_systemSize[2])  atom->position[2] -= m_systemSize[2];

        if (atom->position[0] < 0) atom->position[0] += m_systemSize[0];

        if (atom->position[1] < 0) atom->position[1] += m_systemSize[1];

        if (atom->position[2] < 0) atom->position[2] += m_systemSize[2];
}}

void System::removeTotalMomentum() {
    total_momentum.zeros();
    double numb_atoms = 0;
    for(Atom *atom : m_atoms) {
        total_momentum += atom->velocity * atom->mass();
        numb_atoms++;
    }
    vec3 total_momentum_per_atom = total_momentum / numb_atoms;
    total_momentum.zeros();

    for(Atom *atom : m_atoms) {
        atom->velocity += -1*total_momentum_per_atom / atom->mass();
        total_momentum += atom->velocity * atom->mass();
    }
}
void System::checkMomentum(){
    total_momentum.zeros();
    for(Atom *atom : m_atoms) {
        total_momentum += atom->velocity * atom->mass();
    }
    double tol = 1e-10;
    if (total_momentum.length() > tol ){
        std::cout << "Error: total momentum not conserved" <<std::endl;
        std::cout << "Current net momentum = "<< total_momentum.length() <<std::endl;
        //exit(EXIT_FAILURE);
    }

}
double System::meanDisplacement(){
    meanDispSqr = 0;
    vec3 Displacement;
    for(Atom *atom : m_atoms) {
        Displacement = atom->position - atom->initial_position;
        for (int u = 0; u<3; u++){
            if (Displacement[u] >   0.5*m_systemSize[u]) Displacement[u] += -m_systemSize[u];
            if (Displacement[u] <= -0.5*m_systemSize[u]) Displacement[u] +=  m_systemSize[u];
        }
        meanDispSqr += Displacement.lengthSquared();

    }
    meanDispSqr /= N;
    return meanDispSqr;
}

void System::setPeriodicSize(double PeriodicSize)
{
    Periodicsize = PeriodicSize;
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // Number of atoms
    N = 4*numberOfUnitCellsEachDimension*numberOfUnitCellsEachDimension*numberOfUnitCellsEachDimension;

    double x,y,z;
    for(int i=0; i < numberOfUnitCellsEachDimension; i++){
        for(int j=0; j < numberOfUnitCellsEachDimension; j++){
            for(int k=0; k < numberOfUnitCellsEachDimension; k++){
                Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = i*latticeConstant;
                y = j*latticeConstant;
                z = k*latticeConstant;
                atom->position.set(x,y,z);
                atom->initial_position = atom->position;
                atom->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom);

                Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = i*latticeConstant;
                y = (j+j+1)*latticeConstant/2;
                z = (k+k+1)*latticeConstant/2;
                atom1->position.set(x,y,z);
                atom1->initial_position = atom1->position;
                atom1->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom1);

                Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = (i+i+1)*latticeConstant/2;
                y = j*latticeConstant;
                z = (k+k+1)*latticeConstant/2;
                atom2->position.set(x,y,z);
                atom2->initial_position = atom2->position;
                atom2->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom2);

                Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = (i+i+1)*latticeConstant/2;
                y = (j+j+1)*latticeConstant/2;
                z = k*latticeConstant;
                atom3->position.set(x,y,z);
                atom3->initial_position = atom3->position;
                atom3->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom3);
    }}}
    setSystemSize(vec3((numberOfUnitCellsEachDimension)*latticeConstant,(numberOfUnitCellsEachDimension)*latticeConstant,(numberOfUnitCellsEachDimension)*latticeConstant));
}

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.m_potentialEnergy = 0;
    int i=0, j=0;
    vec3 r_ij_vec;
    double r_ij = 0;
    double sigma_div_r_ij_sqr;
    double sigma_div_r_ij_12;
    double sigma_div_r_ij_6;
    for(Atom *atomi : m_atoms) {
        i++;
        j = 0;
        for(Atom *atomj : m_atoms) {
            j++;
            r_ij_vec = atomj->position - atomi->position;
            for (int u = 0; u<3; u++){
                if (r_ij_vec[u] >   0.5*m_systemSize[u]) r_ij_vec[u] += -m_systemSize[u];
                if (r_ij_vec[u] <= -0.5*m_systemSize[u]) r_ij_vec[u] +=  m_systemSize[u];
            }
            r_ij = r_ij_vec.length();

            sigma_div_r_ij_sqr = (m_potential.sigma() / r_ij)*m_potential.sigma() / r_ij;
            sigma_div_r_ij_12 = (sigma_div_r_ij_sqr)*(sigma_div_r_ij_sqr)*(sigma_div_r_ij_sqr)*(sigma_div_r_ij_sqr)*(sigma_div_r_ij_sqr)*(sigma_div_r_ij_sqr);
            sigma_div_r_ij_6  = (sigma_div_r_ij_sqr)*(sigma_div_r_ij_sqr)*(sigma_div_r_ij_sqr);

            if ( i != j && r_ij != 0) atomi->force += - 24*(m_potential.epsilon()/r_ij)*(2* sigma_div_r_ij_12 - sigma_div_r_ij_6)* (r_ij_vec)/r_ij;

            if ( i > j && i != j && r_ij != 0 ) m_potential.m_potentialEnergy += 4*m_potential.epsilon()*(sigma_div_r_ij_12-sigma_div_r_ij_6);
        }
    }
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);
    m_steps++;
    m_time += dt;
}
