/*
 * Schwarzschild geodesics in the equatorial plane.
 * Integrates a bundle of test particles (massive + photons) with RK4
 * in an affine parameter lambda (= proper time tau for massive).
 * Writes one HDF5 group per orbit for visualization in Python/Manim.
 *
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>

// Units: G = c = 1, mass m sets the length scale (1 M = GM/c^2)

#define M_BH        1.0       // Mass of the black hole
#define dLAMBDA     0.05      // affine-parameter step
#define STRIDE      10        // Save every Nth step
#define R_MAX       60.0      // Stop when r exceeds this (escape)
#define R_HORIZON   (2.0 * M_BH * 1.01)  // Stop just outside horizon (capture)

/* Per-orbit step budgets */
#define N_BOUND     20000     // bound massive — a few orbital periods
#define N_PLUNGE    20000     // plunge — captured well before this
#define N_PHOTON    100000    // photons — ring photons need many steps


/* State vector y = [r, p_r, phi, t]*/
typedef struct {
    double r;      // radial coordinate
    double p_r;    // radial momentum
    double phi;    // angular coordinate
    double t;      // coordinate time
} State;

/* Initial conditions + metadata for one orbit */
typedef struct {
    int epsilon;          // 1 = massive, 0 = photon
    double r0;
    double phi0;
    double p_r0;
    double L;             // conserved angular momentum
    double E;             // conserved energy
    int n_steps;          // integration step budget for this orbit
    const char *kind;     // string tag for Manim coloring
} IC;


/*Right-hand side of the geodesic equations.
 * The dpr_dlambda term: -epsilon*M/r^2 + L^2/r^3 - 3*M*L^2/r^4
*/
static State rhs(State y, int epsilon, double L, double E){
    State dy;
    double r = y.r;
    double f = 1.0 - 2.0*M_BH/r; // Schwarzschild metric factor

    dy.r   = y.p_r;
    dy.p_r = -epsilon*M_BH/pow(r,2) + pow(L,2)/pow(r,3) - 3*M_BH*pow(L,2)/pow(r,4);
    dy.phi = L/pow(r,2);
    dy.t   = E/f;
    return dy;
}

/* RK4 integration step */
static State rk4_step(State y, double h, int epsilon, double L, double E){
    State k1 = rhs(y, epsilon, L, E);
    State y2 = {y.r + 0.5*h*k1.r, y.p_r + 0.5*h*k1.p_r, y.phi + 0.5*h*k1.phi, y.t + 0.5*h*k1.t};

    State k2 = rhs(y2, epsilon, L, E);
    State y3 = {y.r + 0.5*h*k2.r, y.p_r + 0.5*h*k2.p_r, y.phi + 0.5*h*k2.phi, y.t + 0.5*h*k2.t};

    State k3 = rhs(y3, epsilon, L, E);
    State y4 = {y.r + h*k3.r, y.p_r + h*k3.p_r, y.phi + h*k3.phi, y.t + h*k3.t};

    State k4 = rhs(y4, epsilon, L, E);

    State out;
    out.r   = y.r   + (h/6.0)*(k1.r   + 2.0*k2.r   + 2.0*k3.r   + k4.r);
    out.p_r = y.p_r + (h/6.0)*(k1.p_r + 2.0*k2.p_r + 2.0*k3.p_r + k4.p_r);
    out.phi = y.phi + (h/6.0)*(k1.phi + 2.0*k2.phi + 2.0*k3.phi + k4.phi);
    out.t   = y.t   + (h/6.0)*(k1.t   + 2.0*k2.t   + 2.0*k3.t   + k4.t);
    return out;
}

/* Massive particle released at a turning point (p_r = 0).
 * E follows from the constraint: E = sqrt((1-2M/r0)(1 + L^2/r0^2)).
*/
static IC make_massive_bound(double r0, double phi0, double L, const char *kind){
    double f = 1.0 - 2.0*M_BH/r0;
    double E = sqrt(f*(1.0 + pow(L,2)/pow(r0,2)));
    IC ic = {1, r0, phi0, 0.0, L, E, N_BOUND, kind};
    return ic;
}

/* Radial plunge from rest at r0 (L = 0, p_r0 = 0, E = sqrt(1 - 2M/r0)) */
static IC make_radial_plunge(double r0, double phi0, const char *kind){
    double f = 1.0 - 2.0*M_BH/r0;
    double E = sqrt(f);
    IC ic = {1, r0, phi0, 0.0, 0.0, E, N_PLUNGE, kind};
    return ic;
}

/* Photon coming in from the right (positive x) with impact parameter b.
 * Convention: E = 1, L = b, phi0 = asin(b/r0), p_r0 negative (infalling).
*/
static IC make_photon(double r0, double b, const char *kind){
    double phi0 = asin(b/r0);
    double f = 1.0 - 2.0*M_BH/r0;
    double E = 1.0;
    double L = b;
    double V = f * pow(L,2) / pow(r0,2);   // photon V_eff(r0)
    double p_r0 = -sqrt(E*E - V);
    IC ic = {0, r0, phi0, p_r0, L, E, N_PHOTON, kind};
    return ic;
}


/* HDF5 attribute helpers */
static void attr_double(hid_t loc, const char *name, double v){
    hid_t s = H5Screate(H5S_SCALAR);
    hid_t a = H5Acreate(loc, name, H5T_NATIVE_DOUBLE, s, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(a, H5T_NATIVE_DOUBLE, &v);
    H5Aclose(a); H5Sclose(s);
}

static void attr_int(hid_t loc, const char *name, int v){
    hid_t s = H5Screate(H5S_SCALAR);
    hid_t a = H5Acreate(loc, name, H5T_NATIVE_INT, s, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(a, H5T_NATIVE_INT, &v);
    H5Aclose(a); H5Sclose(s);
}

static void attr_string(hid_t loc, const char *name, const char *v){
    hid_t s = H5Screate(H5S_SCALAR);
    hid_t t = H5Tcopy(H5T_C_S1);
    H5Tset_size(t, strlen(v));
    hid_t a = H5Acreate(loc, name, t, s, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(a, t, v);
    H5Aclose(a); H5Tclose(t); H5Sclose(s);
}


/* Integrate one orbit and dump it into its own HDF5 group */
static void integrate_orbit(hid_t file, const char *gname, IC ic){
    State y = {ic.r0, ic.p_r0, ic.phi0, 0.0};

    /* Allocate memory for the output arrays */
    size_t n_save = ic.n_steps / STRIDE + 1;
    double *lambda_arr = malloc(n_save * sizeof(double));
    double *r_arr      = malloc(n_save * sizeof(double));
    double *phi_arr    = malloc(n_save * sizeof(double));
    double *x_arr      = malloc(n_save * sizeof(double)); // Cartesian for plotting
    double *y_arr      = malloc(n_save * sizeof(double));
    double *t_arr      = malloc(n_save * sizeof(double)); // Schwarzschild coordinate time

    /*Integration loop */
    size_t k = 0;
    double lambda = .0;
    lambda_arr[k] = lambda;
    r_arr[k]      = y.r;
    phi_arr[k]    = y.phi;
    x_arr[k]      = y.r * cos(y.phi);
    y_arr[k]      = y.r * sin(y.phi);
    t_arr[k]      = y.t;
    k++;

    int fate = 0; // 0 = ran out of steps, 1 = captured, 2 = escaped
    for (int i = 1; i <= ic.n_steps; i++){
        y = rk4_step(y, dLAMBDA, ic.epsilon, ic.L, ic.E);
        lambda += dLAMBDA;

        if (y.r <= R_HORIZON){ fate = 1; break; }
        if (y.r >= R_MAX)    { fate = 2; break; }

        if (i % STRIDE == 0){
            lambda_arr[k] = lambda;
            r_arr[k]      = y.r;
            phi_arr[k]    = y.phi;
            x_arr[k]      = y.r * cos(y.phi);
            y_arr[k]      = y.r * sin(y.phi);
            t_arr[k]      = y.t;
            k++;
        }
    }

    size_t n_written = k;
    const char *fate_str = (fate == 1) ? "captured" : (fate == 2) ? "escaped" : "timeout";
    printf("[%s] %-8s lambda = %g, n = %zu\n", gname, fate_str, lambda, n_written);

    /* Write group + datasets + attributes */
    hid_t grp = H5Gcreate(file, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dims[1] = {n_written};
    hid_t space = H5Screate_simple(1, dims, NULL);

    #define WRITE(name,buf) do{ \
        hid_t dset = H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); \
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf); \
        H5Dclose(dset); \
    } while(0)

    WRITE("lambda", lambda_arr);
    WRITE("r",      r_arr);
    WRITE("phi",    phi_arr);
    WRITE("x",      x_arr);
    WRITE("y",      y_arr);
    WRITE("t",      t_arr);

    #undef WRITE

    attr_int   (grp, "epsilon", ic.epsilon);
    attr_double(grp, "E",       ic.E);
    attr_double(grp, "L",       ic.L);
    attr_double(grp, "r0",      ic.r0);
    attr_double(grp, "phi0",    ic.phi0);
    attr_double(grp, "p_r0",    ic.p_r0);
    attr_string(grp, "kind",    ic.kind);
    attr_string(grp, "fate",    fate_str);

    H5Sclose(space);
    H5Gclose(grp);

    free(lambda_arr); free(r_arr); free(phi_arr);
    free(x_arr); free(y_arr); free(t_arr);
}

int main(void){
    double b_crit = 3.0 * sqrt(3.0) * M_BH;   // ≈ 5.196152, photon-sphere impact parameter

    /* Orbit table — mix of bound massive, plunges, and photons spanning b_crit */
    IC orbits[] = {
        /* Bound massive — show perihelion precession (distributed around the BH) */
        make_massive_bound(10.0, 0.0, 3.8, "massive_bound"),
        make_massive_bound(12.0, 1.2, 4.0, "massive_bound"),
        make_massive_bound(15.0, 2.5, 4.5, "massive_bound"),
        make_massive_bound( 8.0, 4.0, 3.9, "massive_bound"),
        make_massive_bound(20.0, 5.5, 4.2, "massive_bound"),

        /* Near-ISCO bound orbit — tightly precessing rosette, visually distinct */
        make_massive_bound( 6.5, 2.0, 3.5, "massive_near_isco"),

        /* Radial plunges from rest */
        make_radial_plunge(20.0, 0.7, "massive_plunge"),
        make_radial_plunge(15.0, 3.6, "massive_plunge"),

        /* Photons — captured well away from b_crit */
        make_photon(30.0, 2.0, "photon_capture"),
        make_photon(30.0, 4.0, "photon_capture"),
        make_photon(30.0, 5.0, "photon_capture"),

        /* Photons winding the photon sphere — exponentially close to b_crit */
        make_photon(30.0, b_crit - 1e-2, "photon_capture"),
        make_photon(30.0, b_crit - 1e-3, "photon_capture"),
        make_photon(30.0, b_crit - 1e-5, "photon_capture"),
        make_photon(30.0, b_crit - 1e-7, "photon_capture"),
        make_photon(30.0, b_crit + 1e-7, "photon_escape"),
        make_photon(30.0, b_crit + 1e-5, "photon_escape"),
        make_photon(30.0, b_crit + 1e-3, "photon_escape"),
        make_photon(30.0, b_crit + 1e-2, "photon_escape"),

        /* Photons — escaping with progressively less deflection */
        make_photon(30.0, 5.5,  "photon_escape"),
        make_photon(30.0, 6.0,  "photon_escape"),
        make_photon(30.0, 7.0,  "photon_escape"),
        make_photon(30.0, 8.0,  "photon_escape"),
        make_photon(30.0, 10.0, "photon_escape"),
    };
    int n_orbits = sizeof(orbits) / sizeof(orbits[0]);

    int n_massive = 0, n_photon = 0;
    for (int i = 0; i < n_orbits; i++){
        if (orbits[i].epsilon) n_massive++;
        else                   n_photon++;
    }
    printf("Orbits: %d total (%d massive, %d photons)\n", n_orbits, n_massive, n_photon);

    /* Create file and stamp global attributes on the root group */
    hid_t file = H5Fcreate("geodesic.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0){
        fprintf(stderr, "Failed to create geodesic.h5\n");
        return 1;
    }
    hid_t root = H5Gopen(file, "/", H5P_DEFAULT);
    attr_double(root, "M",         M_BH);
    attr_double(root, "r_horizon", 2.0*M_BH);
    attr_double(root, "r_isco",    6.0*M_BH);
    attr_double(root, "r_photon",  3.0*M_BH);
    attr_double(root, "b_crit",    3.0*sqrt(3.0)*M_BH);
    attr_int   (root, "n_orbits",  n_orbits);
    H5Gclose(root);

    /* Run every orbit into /orbit_NNN */
    for (int i = 0; i < n_orbits; i++){
        char gname[32];
        snprintf(gname, sizeof(gname), "/orbit_%03d", i);
        integrate_orbit(file, gname, orbits[i]);
    }

    H5Fclose(file);

    printf("Wrote %d orbits to geodesic.h5\n", n_orbits);
    return 0;
}
