/**
 * Adding custom post-timestep modifications and forces.
 *
 * This allows the user to use the built-in functions of REBOUNDx
 * but also include their own specialised functions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"


static inline struct reb_particle rebx_particle_minus(struct reb_particle p1, struct reb_particle p2){
    struct reb_particle p = {0};
    p.m = p1.m-p2.m;
    p.x = p1.x-p2.x;
    p.y = p1.y-p2.y;
    p.z = p1.z-p2.z;
    p.vx = p1.vx-p2.vx;
    p.vy = p1.vy-p2.vy;
    p.vz = p1.vz-p2.vz;
    p.ax = p1.ax-p2.ax;
    p.ay = p1.ay-p2.ay;
    p.az = p1.az-p2.az;
    return p;
}
static inline void rebx_subtract_posvel(struct reb_particle* p, struct reb_particle* diff, const double massratio){
    p->x -= massratio*diff->x;
    p->y -= massratio*diff->y;
    p->z -= massratio*diff->z;
    p->vx -= massratio*diff->vx;
    p->vy -= massratio*diff->vy;
    p->vz -= massratio*diff->vz;
}

static struct reb_particle inference_drag_operator(struct reb_simulation* const sim, struct rebx_operator* const dragoperator, const double dt){
    double* dragconst = rebx_get_param(sim->extras, dragoperator->ap, "dragconst");    // get parameters we want user to set
    double* angle = rebx_get_param(sim->extras, dragoperator->ap, "angle");    // get parameters we want user to set
    int err=0;
    struct reb_particle* const p1 = &sim->particles[1];
    struct reb_particle* const primary = &sim->particles[0];
    struct reb_orbit o1 = reb_orbit_from_particle_err(sim->G, sim->particles[1], sim->particles[0], &err);
    struct reb_orbit o2 = reb_orbit_from_particle_err(sim->G, sim->particles[2], sim->particles[0], &err);
    double phi_angle = 0.;
    if(angle != NULL){                              
        phi_angle = (*angle);
    }
    double a1 = o1.a;
    double e1 = o1.e;
    double l1 = o1.l;
    double pom1 = o1.pomega;
    double l2 = o2.l;

    double phi = 2.*l2-l1-pom1+phi_angle;
    
    if(dragconst != NULL){                              
        a1 += 12.*a1*e1*cos(phi)*dt/(*dragconst);
        e1 -= cos(phi)*dt/(*dragconst);
        if (e1<1e-5){
            e1 = 1e-5;
        }
    }
    double m1 = p1->m;
    double GG = sim->G;
    return reb_particle_from_orbit(GG, *primary, m1, a1, e1, o1.inc, o1.Omega, o1.omega, o1.f);
}

void inference_drag(struct reb_simulation* const sim, struct rebx_operator* const dragoperator, const double dt){
    struct reb_particle com = reb_simulation_com(sim); // Start with full com for jacobi and barycentric coordinates.
    struct reb_particle* p = &sim->particles[1];

    struct reb_particle modified_particle = inference_drag_operator(sim, dragoperator, dt);
    struct reb_particle diff = rebx_particle_minus(modified_particle, *p);
    p->x = modified_particle.x;
    p->y = modified_particle.y;
    p->z = modified_particle.z;
    p->vx = modified_particle.vx;
    p->vy = modified_particle.vy;
    p->vz = modified_particle.vz;

    double massratio;
    massratio = p->m/(com.m + p->m);

    for(int j=0; j < 2 + 1; j++){    // stop at j=i if inclusive, at i-1 if not
        rebx_subtract_posvel(&sim->particles[j], &diff, massratio);
    }
}

int run_sim(double tau_nAeA, double tau_nAeB, double tau_nAnB, double mass_ratio, double phi_angle, double tau_nA0){
    char filename[50];
    sprintf(filename, "phi%.2f_%.2f_%.2f_%.2f_%.2f.bin", phi_angle, tau_nAeB, tau_nAnB, mass_ratio, tau_nA0);
    // char* filename = "simulationarchive_inference.bin";
    // parameter definitions:
    
    double tau_nA = 1.e5 *2. *M_PI;
    double tau_eA = tau_nA/tau_nAeA/(1.-1./tau_nAnB);
    double tau_nB = tau_nA*tau_nAnB;
    double tau_eB = tau_nA/tau_nAeB*tau_nAnB/(1.-1./tau_nAnB);
    double tau_0 = -tau_nA/tau_nA0;
    // double tau_0 = INFINITY;
    
    struct reb_simulation* sim = reb_simulation_create();
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->dt = 0.004; 

    struct reb_particle p = {0}; 
    p.m     = 0.31;
    reb_simulation_add(sim, p); 
    reb_simulation_add_fmt(sim, "m a e", 5e-5/(1+mass_ratio), 1.46, 0.0001);
    reb_simulation_add_fmt(sim, "m a e", 5e-5/(1+mass_ratio)*mass_ratio, 2.47, 0.0001);
    reb_simulation_move_to_com(sim); //add particles and move to the center of mass

    struct rebx_extras* rebx = rebx_attach(sim);  // first initialize rebx

    struct rebx_operator* drag = rebx_create_operator(rebx, "inference_drag"); // Now we add inference torque
        drag->operator_type = REBX_OPERATOR_UPDATER;
    drag->step_function = inference_drag;  // set function pointer to what we wrote above
    rebx_add_operator(rebx, drag);      // Now it's initialized, add to REBOUNDx

    struct rebx_force* mo = rebx_load_force(rebx, "modify_orbits_forces"); // add smooth migration & damping
    rebx_add_force(rebx, mo);

    rebx_register_param(rebx, "dragconst", REBX_TYPE_DOUBLE);  // register inference torque parameter
    rebx_set_param_double(rebx, &drag->ap, "dragconst", tau_0); // define the torque magnitude.

    rebx_register_param(rebx, "angle", REBX_TYPE_DOUBLE);  // register inference torque parameter
    rebx_set_param_double(rebx, &drag->ap, "angle", phi_angle); // define the torque magnitude.

    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -tau_nB/2.);         // add semimajor axis damping on inner planet (e-folding timescale)
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_e", -tau_eB); // add linear precession (set precession period). Won't do anything for modify_orbits_forces

    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", -tau_nA/2.);         // add semimajor axis damping on inner planet (e-folding timescale)
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_e", -tau_eA); // add linear precession (set precession period). Won't do anything for modify_orbits_forces


    double tmax = 10.e4/(1.-1./tau_nAnB)/2.5;

    if (remove(filename) == 0) {
        printf("Existing binary file deleted: %s\n", filename);
    } else {
        printf("No existing binary file to delete: %s\n", filename);
    }
    reb_simulation_save_to_file_interval(sim,filename, tmax/5.e3);
    reb_simulation_integrate(sim, tmax);
    rebx_free(rebx);                            // Free all the memory allocated by rebx
}


int main(int argc, char* argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <i> <j> <k> <l>\n", argv[0]);
        return 1;
    }

    double dynamic_tau_nAeB = strtod(argv[1], NULL); //1-100
    double dynamic_param = strtod(argv[2], NULL);
    double mass_ratio = strtod(argv[3], NULL);
    double phi_angle = strtod(argv[4], NULL);
    double tau_nA0 = strtod(argv[5], NULL);

    double tau_nAeA = dynamic_tau_nAeB;
    run_sim(tau_nAeA, dynamic_tau_nAeB, dynamic_param, mass_ratio, phi_angle, tau_nA0);
    return 0;
}