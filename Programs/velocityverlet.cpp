#include "velocityverlet.h"
#include "system.h"
#include "atom.h"
#include "vector"
void VelocityVerlet::integrate(System &system, double dt)
{
    if(m_firstStep) {
        system.calculateForces();
        m_firstStep = false;
    }
    for(Atom *atom : system.atoms()) {
        atom->position += dt*atom->velocity + 0.5*dt*dt*atom->force/atom->mass();
    }
    system.applyPeriodicBoundaryConditions();
    system.calculateForces(); // New positions, recompute forces

    for(Atom *atom : system.atoms()) {
        atom->velocity += 0.5*dt*(atom->force + atom->force_prev)/atom->mass();
        atom->force_prev = atom->force;
    }
}

