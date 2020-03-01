#ifdef _WIN32
#include <direct.h>
#define getcwd _getcwd
#elif linux
#include <unistd.h>
#endif

#define CELL_DIVISION 1.0
#define SKIN 0.0
#define MAX_CELL_NEIG 1000
#define MAX_ATOMS_PER_CELL 500

#define USEOPENMP

#define DONAVIER

#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <omp.h>
//-------------------------------------------------------------------------
using namespace std;

int nthreads;
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//const double kboltz = 0.0019872067;     /* boltzman constant in kcal/mol/K */
const double kboltz = 0.008314462175;     /* boltzman constant in kJ/mol/K */
FILE *logfile, *dumpfile, *pcglog, *navierfile;

struct _navier
{
    int nx, ny, nstep, maxit;
    double dt, mu, beta, dx, dy;
    double **u, **v, **p;
    double **ut, **vt, **c;
    double **uu, **vv, **w;
    double bc_u_north, bc_u_south, bc_v_east, bc_v_west;
    double bc_u_east, bc_u_west;
};
typedef struct _navier navier_t;

/* structure to hold the complete information about the MD system */
struct _mdsys
{
    int natoms, nbonds, nangles, nsteps;
    int ntypes, nbtypes, natypes;
    int nxtout, nveout, step;
    double L, B, H, halfL, halfB, halfH;
    bool error_flag;
    double dt, *mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp, tref;
    int *id, *type, *molid;
    double *charge;
    double *x, *y, *z;
    double *vx, *vy, *vz;
    double *ax, *ay, *az;
    double *fx, *fy, *fz;
    int *bondlist, *anglelist;
    int nthreads;
};
typedef struct _mdsys mdsys_t;

struct _coeff
{
    double lj_A[11][11], lj_B[11][11]; // maximum ten types
    double kbond[11];
    double d0[11];
    double kangle[11];
    double theta0[11];
};
typedef struct _coeff coeff_t;

struct _cellist
{
    int ncells, neighcells;
    int **neighlist;
    int **atomlist;
    int *nlocal;
    double cellsizex, cellsizey, cellsizez;
    int numcellx, numcelly, numcellz;
};
typedef struct _cellist cellist_t;

struct _rcgrelated
{
    int ntype2, ntype3, ntype4;
    int *type2list, *type3list, *type4list;
    //#pair_sumith k1 k2 rbond rcrit
    double k1, rbond, k2, rcrit, rcrit_sqr;
    double pcgbond, pcg_pair_A, pcg_pair_B;
    int **pcgneighlist;
};
typedef struct _rcgrelated rcgmd_t;
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
int READ_DATA_FILE(mdsys_t *, string, navier_t *);
string READ_INPUT_FILE(mdsys_t *, coeff_t *);
int WRITE_DUMP(mdsys_t *, int, rcgmd_t *rcgmd);
void SET_COEFF(int L_lo, int L_hi, int R_lo, int R_hi, double eps, double sigma, coeff_t *coeff);
void velverlet(mdsys_t *, coeff_t *, cellist_t *celist, rcgmd_t *rcgmd);
void force(mdsys_t *, coeff_t *, cellist_t *celist, rcgmd_t *rcgmd);
double pbc(double x, const double );
void ekin(mdsys_t *);
void debug_function(mdsys_t *, int which);
void check_and_update(mdsys_t *);
void initial_momentum(mdsys_t *);
void welcome_message();
void init_rcgmd(mdsys_t *sys, rcgmd_t *rcgmd);
void rcgmd_force(mdsys_t *sys, rcgmd_t *rcgmd);
double heavyside(double val);
void solve_navier(navier_t *);
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
string s_cwd;
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
int main()
{
    string datafile;
    mdsys_t sys;
    coeff_t coeff;
    cellist_t celist;
    rcgmd_t rcgmd;

    navier_t navier;

    time_t tstart, tend;
    tstart = time(0);
    s_cwd = getcwd(NULL,0);
    cout << "CWD is: " << s_cwd << endl;

    sys.error_flag = 0;

    logfile = fopen("logfile.res","w");
    dumpfile = fopen("dump.mol","w");
    pcglog = fopen("pcg_log.log","w");

    welcome_message();

    datafile = READ_INPUT_FILE(&sys, &coeff);
    if( READ_DATA_FILE(&sys, datafile, &navier) == -1)
    {
        fprintf(logfile,"\nError reading data file...exiting\n\n");
        sys.error_flag = 1;
    }

    initial_momentum(&sys);

    check_and_update(&sys);

    init_rcgmd(&sys, &rcgmd);

#ifdef USEOPENMP
    #pragma omp parallel
    {
        if(0 == omp_get_thread_num())
        {
            nthreads=omp_get_num_threads();
            fprintf(logfile,"\n\nRunning RCGMD OpenMP version using %d threads\n\n", nthreads);
        }
    }
#endif // USEOPENMP
    fprintf(logfile,"\nRunning the main program, please wait...\n\n");


    for(sys.step = 1; sys.step <= sys.nsteps; ++sys.step)
    {
        if(sys.error_flag == 0)
        {
            check_and_update(&sys);

#ifdef DONAVIER
            //put a condition here for multi time scale
            solve_navier(&navier);
#endif // DONAVIER

            //debug_function(&sys,1); //debug 1.force, 2.mass, 3.velo, 4disp
            if ((sys.step % sys.nxtout) == 0)
                WRITE_DUMP(&sys, 1, &rcgmd); //traj
            if ((sys.step % sys.nveout) == 0)
            {
                WRITE_DUMP(&sys, 2, &rcgmd); //energy
                WRITE_DUMP(&sys, 3, &rcgmd); //pcglog
            }
            velverlet(&sys, &coeff, &celist, &rcgmd);

            //debug_function(&sys,1); //debug 1.force, 2.mass, 3.velo, 4disp
            //debug_function(&sys,2); //debug 1.force, 2.mass, 3.velo, 4disp
            //debug_function(&sys,3); //debug 1.force, 2.mass, 3.velo, 4disp
            //debug_function(&sys,4); //debug 1.force, 2.mass, 3.velo, 4disp

            ekin(&sys);

            fflush(logfile);
            fflush(dumpfile);
            fflush(pcglog);
        }
        else
        {
            fprintf(logfile, "error: exiting...\n");
            break;
        }
    }

    free(sys.mass);
    free(sys.x);
    free(sys.y);
    free(sys.z);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);
    free(sys.ax);
    free(sys.ay);
    free(sys.az);
    free(sys.charge);
    free(sys.molid);
    free(sys.type);
    free(sys.id);
    free(sys.bondlist);
    free(sys.anglelist);
    //free(celist.atomlist);
    //free(celist.neighlist);
    //free(celist.nlocal);
    free(rcgmd.type2list);
    free(rcgmd.type3list);
    free(rcgmd.type4list);

    free(navier.u);
    free(navier.v);
    free(navier.p);
    free(navier.ut);
    free(navier.vt);
    free(navier.c);
    free(navier.uu);
    free(navier.vv);
    free(navier.w);

    tend = time(0);

    cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;
    cout << "Program completed.";
    fprintf(logfile, "It took %f seconds\n", difftime(tend, tstart));
    fprintf(logfile, "Program completed");
    fflush(logfile);

    fclose(logfile);
    fclose(dumpfile);
    fclose(pcglog);
    return 0;
}
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
void solve_navier(navier_t *navier)
{
    for(int i = 1; i <= (navier->nx + 1); i++ )
    {
        navier->u[i][1] = 2 * navier->bc_u_south - navier->u[i][2];
        navier->u[i][navier->ny+2] = 2 * navier->bc_u_north - navier->u[i][navier->ny+1];
    }

    for(int i = 1; i <= (navier->ny + 1); i++ )
    {
        navier->v[1][i] = 2 * navier->bc_v_west - navier->v[2][i];
        navier->v[navier->nx+2][i] = 2 * navier->bc_v_east - navier->v[navier->nx+1][i];
    }

    for(int i = 2; i <= (navier->ny + 1); i++ )
    {
        navier->u[1][i] = navier->bc_u_west;
        navier->u[navier->nx+1][i] = navier->bc_u_east;
    }

    // -- temporary u velocity
    for(int i = 2; i <= navier->nx; i++ )
    {
        for(int j = 2; j <= navier->ny+1; j++)
        {
            navier->ut[i][j] = navier->u[i][j] + navier->dt * ( -(0.25/navier->dx) * ( pow(( navier->u[i+1][j] + navier->u[i][j] ),2) - pow( ( navier->u[i][j] + \
                               navier->u[i-1][j] ), 2) + ( navier->u[i][j+1] + navier->u[i][j] ) * (navier->v[i+1][j] + \
                                       navier->v[i][j] ) - ( navier->u[i][j] + navier->u[i][j-1] ) * ( navier->v[i+1][j-1] + navier->v[i][j-1] ) )+ \
                               (navier->mu* pow(navier->dx,-2)) * ( navier->u[i+1][j] + navier->u[i-1][j] + navier->u[i][j+1] + navier->u[i][j-1] - 4*navier->u[i][j] ) );

        }
    }

    // -- temporary v velocity
    for(int i = 2; i <= navier->nx+1; i++ )
    {
        for(int j = 2; j <= navier->ny; j++)
        {
            navier->vt[i][j] = navier->v[i][j] + navier->dt*( -(0.25/navier->dy) * ( (navier->u[i][j+1] + navier->u[i][j] ) * ( navier->v[i+1][j]+ \
                               navier->v[i][j] ) - ( navier->u[i-1][j+1] + navier->u[i-1][j] ) * ( navier->v[i][j] + navier->v[i-1][j] )+ \
                               pow(( navier->v[i][j+1] + navier->v[i][j] ),2) - pow(( navier->v[i][j] + navier->v[i][j-1] ),2) ) + \
                               (navier->mu * pow(navier->dy,-2)) * ( navier->v[i+1][j] + navier->v[i-1][j] + navier->v[i][j+1] + navier->v[i][j-1] - 4*navier->v[i][j] ) );
        }
    }

    // -- solve for pressure
    for(int iter = 0; iter < navier->maxit; iter++)
    {
        for(int i = 2; i <= navier->nx+1; i++ )
        {
            for(int j = 2; j <= navier->ny+1; j++)
            {
                navier->p[i][j] = navier->beta * navier->c[i][j] * ( navier->p[i+1][j] + navier->p[i-1][j] + navier->p[i][j+1] + navier->p[i][j-1] - \
                                  (navier->dx/navier->dt) * ( navier->ut[i][j] - navier->ut[i-1][j] + navier->vt[i][j] - navier->vt[i][j-1]))  \
                                  + (1-navier->beta)*navier->p[i][j];
            }
        }
    }

    // -- velocity correction
    for(int i = 2; i <= navier->nx; i++ )
    {
        for(int j = 2; j <= navier->ny+1; j++)
        {
            navier->u[i][j] = navier->ut[i][j] - (navier->dt/navier->dx)*(navier->p[i+1][j] - navier->p[i][j]);
        }
    }
    for(int i = 2; i <= navier->nx+1; i++ )
    {
        for(int j = 2; j <= navier->ny; j++)
        {
            navier->v[i][j]= navier->vt[i][j] - (navier->dt/navier->dy)*(navier->p[i][j+1] - navier->p[i][j]);
        }
    }

    // -- write the results
    for(int i = 1; i <= navier->nx+1; i++ )
    {
        for(int j = 1; j <= navier->ny+1; j++)
        {
            navier->uu[i][j] = 0.5 * (navier->u[i][j+1] + navier->u[i][j]);
            navier->vv[i][j] = 0.5 * (navier->v[i+1][j] + navier->v[i][j]);
            navier->w[i][j] = (navier->u[i][j+1] - navier->u[i][j] - navier->v[i+1][j] + navier->v[i][j] )/ (navier->dx + navier->dy);
        }
    }

    navierfile = fopen("navier_uu.data","w");
    fprintf(navierfile,"%d\t%d\n", navier->nx+1, navier->ny+1);
    for(int i = 1; i <=navier->nx+1; i++)
    {
        for(int j = 1; j <=navier->ny+1; j++)
        {
            fprintf(navierfile,"%f\t", navier->uu[i][j]);
        }
        fprintf(navierfile,"\n" );
    }
    fclose(navierfile);

    navierfile = fopen("navier_vv.data","w");
    fprintf(navierfile,"%d\t%d\n", navier->nx+1, navier->ny+1);
    for(int i = 1; i <=navier->nx+1; i++)
    {
        for(int j = 1; j <=navier->ny+1; j++)
        {
            fprintf(navierfile,"%f\t", navier->vv[i][j]);
        }
        fprintf(navierfile,"\n" );
    }
    fclose(navierfile);

    navierfile = fopen("navier_w.data","w");
    fprintf(navierfile,"%d\t%d\n", navier->nx+1, navier->ny+1);
    for(int i = 1; i <=navier->nx+1; i++)
    {
        for(int j = 1; j <=navier->ny+1; j++)
        {
            fprintf(navierfile,"%f\t", navier->w[i][j]);
        }
        fprintf(navierfile,"\n" );
    }
    fclose(navierfile);
}
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
void force(mdsys_t *sys, coeff_t *coeff, cellist_t *celist, rcgmd_t *rcgmd)
{
    int btype, ith_atom, jth_atom;
    double kbond, d0, rk, r, dr, fbond;
    int i1, i2, i3, atype;
    double k,t0,delx1, dely1, delz1;
    double delx2, dely2, delz2, rsq1, rsq2, r1r2, c, s;
    double dtheta, tk, e_angle, a, a11, a12, a22;
    double f1x, f1y, f1z, f3x, f3y, f3z;

    // clearing the pcg neigh list
    for(int i = 0; i < sys->natoms; i++)
    {
        for(int j = 0; j <=sys->ntypes; j++)
        {
            rcgmd->pcgneighlist[i][j] = 0;
        }
    }
#ifdef USEOPENMP
    //LJ with openMP implementation---------------------------------------
    double epot=0.0;
    int chunksize;//multithread related
    chunksize = ceil(double(sys->natoms)/double(nthreads));
    #pragma omp parallel reduction(+:epot)
    {
        int i, tid, istart, iend, natoms;
        double rsqr,ffac, rc_sqr;
        double irr, ir6, ir12;
        double rx,ry,rz;
        int itype, jtype;
        int imolid, jmolid;
        double *fx, *fy, *fz;

        rc_sqr = sys->rcut * sys->rcut;
        natoms = sys->natoms;
        fx = sys->fx;
        fy = sys->fy;
        fz = sys->fz;
        tid = omp_get_thread_num();

        istart = chunksize*tid;
        iend = chunksize*(tid+1);
        if(iend > natoms ) iend = natoms;

        for( i = istart ; i < iend; ++i)
        {
            for(int j = 0; j < natoms; ++j)
            {
                if(i == j) continue; //same atom skip
                itype  = sys->type[i];
                //imolid = sys->molid[i];
                //jmolid = sys->molid[j];
                //if(imolid == jmolid) continue; // same molecule skip
                // get distance between particle i and j
                rx = pbc(sys->x[i] - sys->x[j], sys->L);
                if(fabs(rx) > sys->rcut) continue;

                ry = pbc(sys->y[i] - sys->y[j], sys->B);
                if(fabs(ry) > sys->rcut) continue;

                rz = pbc(sys->z[i] - sys->z[j], sys->H);
                if(fabs(rz) > sys->rcut) continue;

                rsqr = rx*rx + ry*ry + rz*rz;

                // compute force and energy if within cutoff
                if (rsqr < rc_sqr)
                {
                    jtype  = sys->type[j];

                    irr = 1.0 / rsqr;
                    ir6 = irr * irr * irr;
                    ir12 = ir6 * ir6;
                    if(rsqr <= rcgmd->rcrit_sqr) // for rcgmd functions
                    {
                        rcgmd->pcgneighlist[i][jtype] += 1;
                        //rcgmd->pcgneighlist[j][itype] += 1;
                    }
                    ffac = (12 * coeff->lj_A[itype][jtype] * ir12 - 6 * coeff->lj_B[itype][jtype] * ir6) * irr;
                    epot += 0.5*(coeff->lj_A[itype][jtype] * ir12 - coeff->lj_B[itype][jtype] * ir6);

                    fx[i] += rx * ffac;
                    fy[i] += ry * ffac;
                    fz[i] += rz * ffac;
                    //fprintf(logfile, "DEBUG: %d, %d, %.20f, %.20f, %.20f, %.20f\n", i+1, j+1, ffac, rx, ry, rz);
                }
            }
        }
        #pragma omp barrier
        sys->epot = epot;
        for(i = istart; i < iend; ++i)
        {
            sys->fx[i] = fx[i];
            sys->fy[i] = fy[i];
            sys->fz[i] = fz[i];
        }
    }
#else
    //untouched working LJ without cell list and openmp
    //Lennard jones--------------------------------------------
    double rsqr,ffac, rc_sqr;
    double irr, ir6, ir12;
    double rx,ry,rz;
    int itype, jtype;
    int imolid, jmolid;
    rc_sqr = sys->rcut * sys->rcut;
    sys->epot = 0;

    for(int i = 0 ; i < (sys->natoms) - 1; ++i)
    {
        itype  = sys->type[i];
        //imolid = sys->molid[i];
        for(int j = i+1; j < (sys->natoms); ++j)
        {
            //jmolid = sys->molid[j];
            //if(imolid == jmolid) continue; // same molecule skip
            // get distance between particle i and j
            rx = pbc(sys->x[i] - sys->x[j], sys->L);
            ry = pbc(sys->y[i] - sys->y[j], sys->B);
            rz = pbc(sys->z[i] - sys->z[j], sys->H);
            rsqr = rx*rx + ry*ry + rz*rz;

            // compute force and energy if within cutoff
            if (rsqr < rc_sqr)
            {
                jtype  = sys->type[j];

                irr = 1.0 / rsqr;
                ir6 = irr * irr * irr;
                ir12 = ir6 * ir6;
                if(rsqr <= rcgmd->rcrit_sqr) // for rcgmd functions
                {
                    rcgmd->pcgneighlist[i][jtype] += 1;
                    rcgmd->pcgneighlist[j][itype] += 1;
                }
                ffac = (12 * coeff->lj_A[itype][jtype] * ir12 - 6 * coeff->lj_B[itype][jtype] * ir6) * irr;
                sys->epot += (coeff->lj_A[itype][jtype] * ir12 - coeff->lj_B[itype][jtype] * ir6);

                sys->fx[i] += rx * ffac;
                sys->fy[i] += ry * ffac;
                sys->fz[i] += rz * ffac;

                sys->fx[j] -= rx * ffac;
                sys->fy[j] -= ry * ffac;
                sys->fz[j] -= rz * ffac;
                //fprintf(logfile, "DEBUG: %d, %d, %.20f, %.20f, %.20f, %.20f\n", i+1, j+1, ffac, rx, ry, rz);
            }
            // if(i > 1350) cout<<i<<","<<j<<endl;
        }
    }
#endif // USEOPENMP

#ifdef USEOPENMP
    double rx,ry,rz,rsqr;
#endif // USEOPENMP

    //Bond forces--------------------------------------------
    for(int i = 0; i < sys->nbonds; i++)
    {
        btype = 1;
        kbond = coeff->kbond[btype];
        d0 = coeff->d0[btype];

        ith_atom = sys->bondlist[i];
        jth_atom = sys->bondlist[i+sys->nbonds];

        rx = pbc(sys->x[ith_atom] - sys->x[jth_atom], sys->L);
        ry = pbc(sys->y[ith_atom] - sys->y[jth_atom], sys->B);
        rz = pbc(sys->z[ith_atom] - sys->z[jth_atom], sys->H);

        rsqr = rx*rx + ry*ry + rz*rz;

        r = sqrt(rsqr);
        dr = r - d0;
        rk = kbond * dr;
        fbond = -2.0 * rk / r;
        sys->epot += rk * dr;

        sys->fx[ith_atom] += rx * fbond;
        sys->fy[ith_atom] += ry * fbond;
        sys->fz[ith_atom] += rz * fbond;

        sys->fx[jth_atom] -= rx * fbond;
        sys->fy[jth_atom] -= ry * fbond;
        sys->fz[jth_atom] -= rz * fbond;
    }

    //Angle forces-----------------------------------------------
    for(int i = 0; i < sys->nangles; i++)
    {
        i1 = sys->anglelist[i];
        i2 = sys->anglelist[i + sys->nangles];
        i3 = sys->anglelist[i + sys->nangles + sys->nangles];

        atype = 1;
        k = coeff->kangle[atype];
        t0 =coeff->theta0[atype];

        // 1st bond
        delx1 = pbc(sys->x[i1] - sys->x[i2], sys->L);
        dely1 = pbc(sys->y[i1] - sys->y[i2], sys->B);
        delz1 = pbc(sys->z[i1] - sys->z[i2], sys->H);
        rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;

        // 1st bond
        delx2 = pbc(sys->x[i3] - sys->x[i2], sys->L);
        dely2 = pbc(sys->y[i3] - sys->y[i2], sys->B);
        delz2 = pbc(sys->z[i3] - sys->z[i2], sys->H);
        rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;

        r1r2 = sqrt(rsq1 * rsq2);

        //# angle (cos and sin)
        c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
        c /= r1r2;

        if(c > 1.0) c = 1.0;
        if(c < -1.0) c = -1.0;

        s = sqrt(1.0 - c * c);
        if (s < 0.001) s = 0.001;
        s = 1.0 / s;

        //# ------ force & energy
        dtheta = acos(c) - t0;
        tk = k * dtheta;
        e_angle = e_angle + tk * dtheta;

        a = -2.0 * tk * s;
        a11 = a * c / rsq1;
        a12 = -a / r1r2;
        a22 = a * c / rsq2;

        f1x = a11 * delx1 + a12 * delx2;
        f1y = a11 * dely1 + a12 * dely2;
        f1z = a11 * delz1 + a12 * delz2;
        f3x = a22 * delx2 + a12 * delx1;
        f3y = a22 * dely2 + a12 * dely1;
        f3z = a22 * delz2 + a12 * delz1;

        //# apply force to each of 3 atoms
        sys->fx[i1] += f1x;
        sys->fy[i1] += f1y;
        sys->fz[i1] += f1z;

        sys->fx[i2] -= f1x + f3x;
        sys->fy[i2] -= f1y + f3y;
        sys->fz[i2] -= f1z + f3z;

        sys->fx[i3] += f3x;
        sys->fy[i3] += f3y;
        sys->fz[i3] += f3z;
    }
}
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
/* compute kinetic energy */
void ekin(mdsys_t *sys)
{
    int i;
    double mass, scale;

    sys->ekin=0.0;

    for (i=0; i < sys->natoms; ++i)
    {
        mass = sys->mass[sys->type[i]];
        sys->ekin += mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->ekin = sys->ekin / 2.0;
    sys->temp = 2.0 * sys->ekin /(3.0 * sys->natoms * kboltz);

    //scale the temperature
    scale = sqrt(sys->tref / sys->temp);
    for (i=0; i < sys->natoms; ++i)
    {
        sys->vx[i] = sys->vx[i] * scale;
        sys->vy[i] = sys->vy[i] * scale;
        sys->vz[i] = sys->vz[i] * scale;
    }
}
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
/* velocity verlet */
void velverlet(mdsys_t *sys, coeff_t *coeff, cellist_t *celist, rcgmd_t *rcgmd)
{
    int i;
    double mass;
    double half_dt;
    half_dt = 0.5 * sys->dt;

    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i < sys->natoms; ++i)
    {
        mass = sys->mass[sys->type[i]];

        sys->x[i] += sys->dt * sys->vx[i];
        sys->y[i] += sys->dt * sys->vy[i];
        sys->z[i] += sys->dt * sys->vz[i];

        sys->vx[i] += half_dt * sys->ax[i];
        sys->vy[i] += half_dt * sys->ay[i];
        sys->vz[i] += half_dt * sys->az[i];

        sys->fx[i] = 0;
        sys->fy[i] = 0;
        sys->fz[i] = 0;
    }
    /* compute forces and potential energy */
    force(sys, coeff, celist, rcgmd);

    rcgmd_force(sys, rcgmd); //this should call after force only

    /* second part: propagate velocities by another half step */
    for (i=0; i<sys->natoms; ++i)
    {
        mass = sys->mass[sys->type[i]];

        sys->ax[i] = sys->fx[i] / mass;
        sys->ay[i] = sys->fy[i] / mass;
        sys->az[i] = sys->fz[i] / mass;

        sys->vx[i] += half_dt * sys->ax[i];
        sys->vy[i] += half_dt * sys->ay[i];
        sys->vz[i] += half_dt * sys->az[i];
    }
}

/* helper function: apply minimum image convention */
double pbc(double x, const double box)
{
    if(fabs(x) > box/2.0)
    {
        if (x < 0)
            x = x + box;
        else
            x = x - box;
    }
    return x;
}

void check_and_update(mdsys_t *sys)
{
    for (int i = 0; i < sys->natoms; ++i)
    {
        if(sys->x[i] > sys->halfL) sys->x[i] = sys->x[i] - sys->L;
        if(sys->y[i] > sys->halfB) sys->y[i] = sys->y[i] - sys->B;
        if(sys->z[i] > sys->halfH) sys->z[i] = sys->z[i] - sys->H;

        if(sys->x[i] < -sys->halfL) sys->x[i] = sys->x[i] + sys->L;
        if(sys->y[i] < -sys->halfB) sys->y[i] = sys->y[i] + sys->B;
        if(sys->z[i] < -sys->halfH) sys->z[i] = sys->z[i] + sys->H;
    }
}
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
int READ_DATA_FILE(mdsys_t *sys, string datafile, navier_t *navier )
{
    string linebuffer, filename, delim;
    filename = "lammps.data";

    ifstream in;
    int length, pos, pos1, pos2;

    in.open(filename.c_str());
    if(!in)
    {
        cout << "Error: unable to open file " <<filename<< endl;
        sys->error_flag = 1;
        return -1;
    }
    else cout << "Reading file " <<filename <<"  pass 1..."<< endl;


    sys->natypes = 0;
    sys->nbtypes = 0;
    sys->nangles = 0;
    sys->nbonds = 0;
    //pass 1 gather data
    while ( getline(in, linebuffer))
    {
        if (linebuffer.find ("atoms") != string::npos)
        {
            length = linebuffer.length();
            pos = linebuffer.find("atoms");
            sys->natoms = atof(linebuffer.substr(0,pos).c_str());
            cout  << "natoms : "<< sys->natoms << endl;
        }
        if (linebuffer.find ("atom types") != string::npos)
        {
            length = linebuffer.length();
            pos = linebuffer.find("atom types");
            sys->ntypes = atof(linebuffer.substr(0,pos).c_str());
            cout  << "ntypes : "<< sys->ntypes << endl;
        }
        if (linebuffer.find ("bonds") != string::npos)
        {
            length = linebuffer.length();
            pos = linebuffer.find("bonds");
            sys->nbonds = atof(linebuffer.substr(0,pos).c_str());
            cout  << "nbonds : "<< sys->nbonds << endl;
        }
        if (linebuffer.find ("bond types") != string::npos)
        {
            length = linebuffer.length();
            pos = linebuffer.find("bond types");
            sys->nbtypes = atof(linebuffer.substr(0,pos).c_str());
            cout  << "nbtypes : "<< sys->nbtypes << endl;
        }
        if (linebuffer.find ("angles") != string::npos)
        {
            length = linebuffer.length();
            pos = linebuffer.find("angles");
            sys->nangles = atof(linebuffer.substr(0,pos).c_str());
            cout  << "nangles : "<< sys->nangles << endl;
        }
        if (linebuffer.find ("angle types") != string::npos)
        {
            length = linebuffer.length();
            pos = linebuffer.find("angle types");
            sys->natypes = atof(linebuffer.substr(0,pos).c_str());
            cout  << "natypes : "<< sys->natypes << endl;
        }
        if (linebuffer.find ("xlo") != string::npos)
        {
            length = linebuffer.length();
            pos2 = linebuffer.find("xlo");
            pos1 = linebuffer.find(" ");
            sys->L = atof(linebuffer.substr(pos1+1, pos2-pos1).c_str()) - atof(linebuffer.substr(0,pos1).c_str());
            sys->halfL = sys->L/2;
            cout  << "L : "<< sys->L << endl;
        }
        if (linebuffer.find ("ylo") != string::npos)
        {
            length = linebuffer.length();
            pos2 = linebuffer.find("ylo");
            pos1 = linebuffer.find(" ");
            sys->B = atof(linebuffer.substr(pos1+1, pos2-pos1).c_str()) - atof(linebuffer.substr(0,pos1).c_str());
            sys->halfB = sys->B/2;
            cout  << "B : "<< sys->B << endl;
        }
        if (linebuffer.find ("zlo") != string::npos)
        {
            length = linebuffer.length();
            pos2 = linebuffer.find("zlo");
            pos1 = linebuffer.find(" ");
            sys->H = atof(linebuffer.substr(pos1+1, pos2-pos1).c_str()) - atof(linebuffer.substr(0,pos1).c_str());
            sys->halfH = sys->H/2;
            cout  << "H : "<< sys->H << endl;
        }
    }
    in.close();

    //allocate memory
    sys->mass=(double *)malloc((sys->ntypes+1)*sizeof(double));
    sys->x = (double *)malloc(sys->natoms *sizeof(double));
    sys->y = (double *)malloc(sys->natoms *sizeof(double));
    sys->z = (double *)malloc(sys->natoms *sizeof(double));
    sys->vx = (double *)malloc(sys->natoms *sizeof(double));
    sys->vy = (double *)malloc(sys->natoms *sizeof(double));
    sys->vz = (double *)malloc(sys->natoms *sizeof(double));
    sys->fx = (double *)malloc(sys->natoms *sizeof(double));
    sys->fy = (double *)malloc(sys->natoms *sizeof(double));
    sys->fz = (double *)malloc(sys->natoms *sizeof(double));
    sys->ax = (double *)malloc(sys->natoms *sizeof(double));
    sys->ay = (double *)malloc(sys->natoms *sizeof(double));
    sys->az = (double *)malloc(sys->natoms *sizeof(double));

    sys->charge = (double *)malloc(sys->natoms *sizeof(double));
    sys->molid = (int *)malloc(sys->natoms *sizeof(int));
    sys->type = (int *)malloc(sys->natoms *sizeof(int));
    sys->id = (int *)malloc(sys->natoms *sizeof(int));
    if(sys->nbonds >0) sys->bondlist = (int *)malloc(sys->nbonds*2 *sizeof(int));
    if(sys->nangles >0) sys->anglelist = (int *)malloc(sys->nangles*3 *sizeof(int));


    for(int i = 0; i < sys->natoms; i++)
    {
        sys->x[i] = sys->y[i] = sys->z[i] = 0;
        sys->vx[i] = sys->vy[i] = sys->vz[i] = 0;
        sys->ax[i] = sys->ay[i] = sys->az[i] = 0;
        sys->fx[i] = sys->fy[i] = sys->fz[i] = 0;
        sys->charge[i] = 0;
    }

    in.open(filename.c_str());
    cout << "Reading mass from " <<filename <<"  pass 2..."<< endl;
    int counter = 0;
    bool started = false;
    //pass 2 gather mass
    while ( getline(in, linebuffer))
    {
        if (linebuffer.find ("Masses") != string::npos)
        {
            started = true;
            continue;
        }
        if (started == true && linebuffer.find (" ") != string::npos && linebuffer.length() > 4)
        {
            counter++;
            length = linebuffer.length();
            pos = linebuffer.find(" ");
            sys->mass[counter] = atof(linebuffer.substr(pos+1, length-pos-1).c_str());
            cout  << "mass : "<< sys->mass[counter] << endl;
        }
        if(counter == sys->ntypes) break;
    }
    in.close();


    in.open(filename.c_str());
    cout << "Reading atoms from " <<filename <<"  pass 3..."<< endl;
    counter = 0;
    started = false;
    //pass 3 gather atom positions
    while ( getline(in, linebuffer))
    {
        if (linebuffer.find ("Atoms") != string::npos)
        {
            started = true;
            continue;
        }
        if (started == true && linebuffer.length() > 2)
        {
            if (linebuffer.find (" ") != string::npos)
                delim = " ";
            else if (linebuffer.find ("\t") != string::npos)
                delim = "\t";
            else
            {
                sys->error_flag = true;
                return -1;
            }

            length = linebuffer.length();
            pos = linebuffer.find(delim);
            sys->id[counter] = atof(linebuffer.substr(0, pos).c_str());
            linebuffer.erase(0, pos+1);

            pos = linebuffer.find(delim);
            sys->molid[counter] = atof(linebuffer.substr(0, pos).c_str());
            linebuffer.erase(0, pos+1);

            pos = linebuffer.find(delim);
            sys->type[counter] = atof(linebuffer.substr(0, pos).c_str());
            linebuffer.erase(0, pos+1);

            pos = linebuffer.find(delim);
            sys->charge[counter] = atof(linebuffer.substr(0, pos).c_str());
            linebuffer.erase(0, pos+1);

            pos = linebuffer.find(delim);
            sys->x[counter] = atof(linebuffer.substr(0, pos).c_str());
            linebuffer.erase(0, pos+1);

            pos = linebuffer.find(delim);
            sys->y[counter] = atof(linebuffer.substr(0, pos).c_str());
            linebuffer.erase(0, pos+1);

            length = linebuffer.length();
            sys->z[counter] = atof(linebuffer.substr(0,length).c_str());

            //fprintf(logfile,"%d, %d, %d, %f, %f, %f, %f\n", sys->id[counter],sys->molid[counter],sys->type[counter],sys->charge[counter],sys->x[counter],sys->y[counter],sys->z[counter] );

            counter++;
        }
        if(counter >= sys->natoms) break;
    }
    in.close();



    if(sys->nbonds > 0)
    {
        in.open(filename.c_str());
        cout << "Reading bonds from " <<filename <<"  pass 4..."<< endl;
        counter = 0;
        started = false;
        //pass 4 gather bonds
        while ( getline(in, linebuffer))
        {
            if (linebuffer.find ("Bonds") != string::npos)
            {
                started = true;
                continue;
            }
            if (started == true && linebuffer.length() > 2)
            {
                if (linebuffer.find (" ") != string::npos)
                    delim = " ";
                else if (linebuffer.find ("\t") != string::npos)
                    delim = "\t";
                else
                {
                    sys->error_flag = true;
                    return -1;
                }

                pos = linebuffer.find(delim);
                linebuffer.erase(0, pos+1);
                pos = linebuffer.find(delim);
                linebuffer.erase(0, pos+1);

                pos = linebuffer.find(delim);
                sys->bondlist[counter] = atof(linebuffer.substr(0, pos).c_str())-1;
                linebuffer.erase(0, pos+1);
                length = linebuffer.length();
                sys->bondlist[counter + sys->nbonds] = atof(linebuffer.substr(0,length).c_str())-1;

                //fprintf(logfile,"%d, %d\n", sys->bondlist[counter],sys->bondlist[counter + sys->nbonds]);

                counter++;
            }
            if(counter >= sys->nbonds) break;
        }
        in.close();
    }

    if(sys->nangles > 0)
    {
        in.open(filename.c_str());
        cout << "Reading angles from " <<filename <<"  pass 5..."<< endl;
        counter = 0;
        started = false;
        //pass 5 gather angles
        while ( getline(in, linebuffer))
        {
            if (linebuffer.find ("Angles") != string::npos)
            {
                started = true;
                continue;
            }
            if (started == true && linebuffer.length() > 2)
            {
                if (linebuffer.find (" ") != string::npos)
                    delim = " ";
                else if (linebuffer.find ("\t") != string::npos)
                    delim = "\t";
                else
                {
                    sys->error_flag = true;
                    return -1;
                }

                pos = linebuffer.find(delim);
                linebuffer.erase(0, pos+1);
                pos = linebuffer.find(delim);
                linebuffer.erase(0, pos+1);

                pos = linebuffer.find(delim);
                sys->anglelist[counter] = atof(linebuffer.substr(0, pos).c_str())-1;
                linebuffer.erase(0, pos+1);

                pos = linebuffer.find(delim);
                sys->anglelist[counter+ sys->nangles] = atof(linebuffer.substr(0, pos).c_str())-1;
                linebuffer.erase(0, pos+1);

                length = linebuffer.length();
                sys->anglelist[counter + 2*sys->nangles] = atof(linebuffer.substr(0,length).c_str())-1;

                //fprintf(logfile,"%d, %d, %d\n", sys->anglelist[counter],sys->anglelist[counter + sys->nangles], sys->anglelist[counter + 2*sys->nangles]);

                counter++;
            }
            if(counter >= sys->nangles) break;
        }
        in.close();
    }

    //---Navier related stuff
    navier->nx = 9;
    navier->ny = 9;
    navier->dt = 0.02;
    navier->nstep = 1;
    navier->mu = 0.1;
    navier->maxit = 100;
    navier->beta = 1.2;
    navier->dx = 1.0 / double(navier->nx); // here 1 is the length
    navier->dy = 1.0 / double(navier->ny); // here 1 is the breadth

    double **dummy = (double **) malloc((navier->nx + 3) *sizeof(double*));
    for(int i = 0; i < (navier->nx + 3); i++)
        dummy[i] = (double *) malloc((navier->ny + 3) * sizeof(double));
    navier->u = dummy;

    double **dummy1 = (double **) malloc((navier->nx + 3) *sizeof(double*));
    for(int i = 0; i <  (navier->nx + 3); i++)
        dummy1[i] = (double *) malloc((navier->ny + 3) * sizeof(double));
    navier->v = dummy1;

    double **dummy2 = (double **) malloc((navier->nx + 3) *sizeof(double*));
    for(int i = 0; i <  (navier->nx + 3); i++)
        dummy2[i] = (double *) malloc((navier->ny + 3) * sizeof(double));
    navier->p = dummy2;

    double **dummy3 = (double **) malloc((navier->nx + 3) *sizeof(double*));
    for(int i = 0; i <  (navier->nx + 3); i++)
        dummy3[i] = (double *) malloc((navier->ny + 3) * sizeof(double));
    navier->ut = dummy3;

    double **dummy4 = (double **) malloc((navier->nx + 3) *sizeof(double*));
    for(int i = 0; i <  (navier->nx + 3); i++)
        dummy4[i] = (double *) malloc((navier->ny + 3) * sizeof(double));
    navier->vt = dummy4;

    double **dummy5 = (double **) malloc((navier->nx + 3) *sizeof(double*));
    for(int i = 0; i <  (navier->nx + 3); i++)
        dummy5[i] = (double *) malloc((navier->ny + 3) * sizeof(double));
    navier->c = dummy5;

    double **dummy6 = (double **) malloc((navier->nx + 3) *sizeof(double*));
    for(int i = 0; i <  (navier->nx + 3); i++)
        dummy6[i] = (double *) malloc((navier->ny + 3) * sizeof(double));
    navier->uu = dummy6;

    double **dummy7 = (double **) malloc((navier->nx + 3) *sizeof(double*));
    for(int i = 0; i <  (navier->nx + 3); i++)
        dummy7[i] = (double *) malloc((navier->ny + 3) * sizeof(double));
    navier->vv = dummy7;

    double **dummy8 = (double **) malloc((navier->nx + 3) *sizeof(double*));
    for(int i = 0; i <  (navier->nx + 3); i++)
        dummy8[i] = (double *) malloc((navier->ny + 3) * sizeof(double));
    navier->w = dummy8;

    for(int i = 0; i < navier->nx+3; i++)
    {
        for(int j = 0; j < navier->ny+3; j++)
        {
            navier->u[i][j] = 0;
            navier->v[i][j] = 0;
            navier->p[i][j] = 0;

            navier->ut[i][j] = 0;
            navier->vt[i][j] = 0;
            navier->c[i][j] = 0.25; //notice the change

            navier->uu[i][j] = 0;
            navier->vv[i][j] = 0;
            navier->w[i][j] = 0;
        }
    }

    //-----------set the conditions for pressure---
    for(int i = 3; i <= navier->ny; i++)
    {
        navier->c[2][i] = 1.0/3.0;
        navier->c[navier->nx + 1][i] = 1.0/3.0;
    }

    for(int i = 3; i <= navier->ny; i++)
    {
        navier->c[i][2] = 1.0/3.0;
        navier->c[i][navier->ny + 1] = 1.0/3.0;
    }

    navier->c[2][2] = 1.0/2.0;
    navier->c[2][navier->ny + 1] = 1.0/2.0;
    navier->c[navier->nx + 1][2] = 1.0/2.0;
    navier->c[navier->nx + 1][navier->ny + 1] = 1.0/2.0;

    //-----------set the boundary conditions for velocities---
    navier->bc_u_north = 0;
    navier->bc_u_south = 0;
    navier->bc_v_east = 0;
    navier->bc_v_west = 0;

    navier->bc_u_east = 0.1; //right
    navier->bc_u_west = 0.1; //left

    return 0;
}
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
string READ_INPUT_FILE(mdsys_t *sys, coeff_t *coeff)
{
    string linebuffer, filename, buffer, datafile_temp, logfile_temp;
    int L_lo, L_hi, R_lo, R_hi, pos, length, bondtype;
    double eps, sigma, kbond, d0;

#ifdef _WIN32
    filename = "conf.in";
#elif linux
    filename= s_cwd + "/conf.in";
#endif

    ifstream in;
    in.open(filename.c_str());
    fprintf(logfile, "Reading input file..\n\n");

    //pass 5 gather angles
    while ( getline(in, linebuffer))
    {
        if(linebuffer.find ("#") != string::npos) continue; // comment line skip it.

        if (linebuffer.find ("log") != string::npos)
        {
            pos = linebuffer.find(" ");
            length = linebuffer.length();
            logfile_temp = linebuffer.substr(pos+1, length - pos).c_str();
            continue;
        }

        if (linebuffer.find ("pair_coeff") != string::npos)
        {
            pos = linebuffer.find(" ");
            buffer = linebuffer.substr(pos+1, 1).c_str();
            if(buffer.compare("*") == 0)
            {
                L_lo = 1;
                L_hi = 10;
            }
            else
            {
                L_lo = L_hi = atof(linebuffer.substr(pos+1, 1).c_str());
            }
            linebuffer.erase(0, pos+1);

            pos = linebuffer.find(" ");
            buffer = linebuffer.substr(pos+1, 1).c_str();
            if(buffer.compare("*") == 0)
            {
                R_lo = 1;
                R_hi = 10;
            }
            else
            {
                R_lo = R_hi = atof(linebuffer.substr(pos+1, 1).c_str());
            }
            linebuffer.erase(0, pos+1);

            pos = linebuffer.find(" ");
            length = linebuffer.length();
            eps = atof(linebuffer.substr(pos+1, length-pos).c_str());
            linebuffer.erase(0, pos+1);

            pos = linebuffer.find(" ");
            length = linebuffer.length();
            sigma = atof(linebuffer.substr(pos+1, length-pos).c_str());
            linebuffer.erase(0, pos+1);

            fprintf(logfile,"pair: %d, %d, %d, %d, %f, %f\n", L_lo, L_hi, R_lo, R_hi,eps,sigma);
            SET_COEFF(L_lo, L_hi, R_lo, R_hi,eps,sigma, coeff);
            continue;
        }

        if (linebuffer.find ("bond_coeff") != string::npos)
        {
            pos = linebuffer.find(" ");
            linebuffer.erase(0, pos+1); //stripped "bond_coeff"
            pos = linebuffer.find(" ");
            bondtype = atof(linebuffer.substr(0, pos).c_str());
            linebuffer.erase(0, pos+1); //stripped bondtype
            pos = linebuffer.find(" ");
            length = linebuffer.length();

            kbond = atof(linebuffer.substr(0, pos).c_str());
            d0 = atof(linebuffer.substr(pos+1, length - pos).c_str());

            coeff->kbond[bondtype] = kbond;
            coeff->d0[bondtype] = d0;
            continue;
        }

        if (linebuffer.find ("angle_coeff") != string::npos)
        {
            pos = linebuffer.find(" ");
            linebuffer.erase(0, pos+1); //stripped "bond_coeff"
            pos = linebuffer.find(" ");
            bondtype = atof(linebuffer.substr(0, pos).c_str());
            linebuffer.erase(0, pos+1); //stripped a
            pos = linebuffer.find(" ");
            length = linebuffer.length();

            kbond = atof(linebuffer.substr(0, pos).c_str());
            d0 = atof(linebuffer.substr(pos+1, length - pos).c_str());

            coeff->kangle[bondtype] = kbond;
            coeff->theta0[bondtype] = d0 * 0.01745329251994329576923690768489; //converting deg to rad
            continue;
        }

        if (linebuffer.find ("timestep") != string::npos)
        {
            pos = linebuffer.find(" ");
            length = linebuffer.length();
            sys->dt = atof(linebuffer.substr(pos+1, length - pos).c_str());
            continue;
        }
        if (linebuffer.find ("rcut") != string::npos)
        {
            pos = linebuffer.find(" ");
            length = linebuffer.length();
            sys->rcut = atof(linebuffer.substr(pos+1, length - pos).c_str());
            continue;
        }
        if (linebuffer.find ("tref") != string::npos)
        {
            pos = linebuffer.find(" ");
            length = linebuffer.length();
            sys->tref = atof(linebuffer.substr(pos+1, length - pos).c_str());
            continue;
        }
        if (linebuffer.find ("nxtout") != string::npos)
        {
            pos = linebuffer.find(" ");
            length = linebuffer.length();
            sys->nxtout = atof(linebuffer.substr(pos+1, length - pos).c_str());
            continue;
        }

        if (linebuffer.find ("nveout") != string::npos)
        {
            pos = linebuffer.find(" ");
            length = linebuffer.length();
            sys->nveout = atof(linebuffer.substr(pos+1, length - pos).c_str());
            continue;
        }

        if (linebuffer.find ("run") != string::npos)
        {
            pos = linebuffer.find(" ");
            length = linebuffer.length();
            sys->nsteps = atof(linebuffer.substr(pos+1, length - pos).c_str());
            continue;
        }
    }
    in.close();
    return datafile_temp;
}

void SET_COEFF(int L_lo, int L_hi, int R_lo, int R_hi, double eps, double sigma, coeff_t *coeff)
{
    for(int i = L_lo; i <= L_hi; ++i)
    {
        for(int j = R_lo; j <= R_hi; j++)
        {
            coeff->lj_A[i][j] = 4 * eps * pow(sigma, 12);
            coeff->lj_B[i][j] = 4 * eps * pow(sigma, 6);

            coeff->lj_A[j][i] = 4 * eps * pow(sigma, 12);
            coeff->lj_B[j][i] = 4 * eps * pow(sigma, 6);
        }
    }
}
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
int WRITE_DUMP(mdsys_t *sys, int which, rcgmd_t *rcgmd)
{
    if(which == 1) //traj
    {
        fprintf(dumpfile,"ITEM: TIMESTEP\n");
        fprintf(dumpfile,"%d\n", sys->step);
        fprintf(dumpfile,"ITEM: NUMBER OF ATOMS\n");
        fprintf(dumpfile,"%d\n", sys->natoms);
        fprintf(dumpfile,"ITEM: BOX BOUNDS pp pp pp\n");
        fprintf(dumpfile,"%f %f\n", -sys->halfL, sys->halfL);
        fprintf(dumpfile,"%f %f\n", -sys->halfB, sys->halfB);
        fprintf(dumpfile,"%f %f\n", -sys->halfH, sys->halfH);
        fprintf(dumpfile,"ITEM: ATOMS id mol type x y z\n");
        for(int i = 0; i < sys->natoms; i++)
        {
            fprintf(dumpfile,"%d %d %d %.3f %.3f %.3f\n", i+1, sys->molid[i], sys->type[i], sys->x[i], sys->y[i], sys->z[i]);
        }
    }
    else if(which == 2) //energy
    {
        if(sys->step == sys->nxtout)
        {
            fprintf(logfile,"step\tke\tpe\ttemp\n");
        }
        fprintf(logfile,"%d\t%.3f\t%.3f\t%.3f\n", sys->step, sys->ekin, sys->epot, sys->temp);
    }
    else if(which == 3) //pcglog
    {
        int pcg_n_E_4 = 0;
        int pcg_n_D_4 = 0;
        int iatom = -1;
        if(sys->step == sys->nxtout)
        {
            fprintf(pcglog,"step\tnE4\tnD4\n");
        }
        for(int i = 0; i < rcgmd->ntype4 ; i++)
        {
            iatom = rcgmd->type4list[i];
            pcg_n_E_4 += rcgmd->pcgneighlist[iatom][3];
            pcg_n_D_4 += rcgmd->pcgneighlist[iatom][2];
        }
        fprintf(pcglog,"%d\t%d\t%d\n", sys->step, pcg_n_E_4, pcg_n_D_4);
    }
    else
    {
        fprintf(logfile, "Error: unknown dump argument");
        sys->error_flag =1;
    }
    return 0;
}
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
void debug_function(mdsys_t *sys, int which)
{
    //debug force values

    if(which== 1)
    {
        fprintf(logfile, "Debugging Force values\n");
        fprintf(logfile, "atomid, FX, FY, FZ\n");
        for(int i = 0; i < sys->natoms; i++)
        {
            fprintf(logfile, "DEBUG: %d, %.20f, %.20f, %.20f\n", i+1, sys->fx[i], sys->fy[i], sys->fz[i]);
        }
        fprintf(logfile, "Debugging complete\n");
    }
    //debug mass values
    if(which== 2)
    {
        fprintf(logfile, "Debugging mass values\n");
        fprintf(logfile, "atomid, type, mass\n");
        for(int i = 0; i < sys->natoms; i++)
        {
            fprintf(logfile, "DEBUG: %d, %d, %f\n", i+1, sys->type[i], sys->mass[sys->type[i]]);
        }
        fprintf(logfile, "Debugging complete\n");
    }
    //debug velocity values
    if(which== 3)
    {
        fprintf(logfile, "Debugging Velocity values\n");
        fprintf(logfile, "atomid, VX, VY, VZ\n");
        for(int i = 0; i < sys->natoms; i++)
        {
            fprintf(logfile, "DEBUG: %d, %f, %f, %f\n", i+1, sys->vx[i], sys->vy[i], sys->vz[i]);
        }
        fprintf(logfile, "Debugging complete\n");
    }
    //debug displacement values
    if(which== 4)
    {
        fprintf(logfile, "Debugging displacement values\n");
        fprintf(logfile, "atomid, x, y, z\n");
        for(int i = 0; i < sys->natoms; i++)
        {
            fprintf(logfile, "DEBUG: %d, %f, %f, %f\n", i+1, sys->x[i], sys->y[i], sys->z[i]);
        }
        fprintf(logfile, "Debugging complete\n");
        fflush(logfile);
    }
    fflush(logfile);
}
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
void initial_momentum(mdsys_t *sys)
{
    srand(time(0));
    double t1, t2, t3;

    for(int i = 0; i < sys->natoms ; i++)
    {
        t1 = rand() / double(RAND_MAX);
        t2 = rand() / double(RAND_MAX);
        t3 = rand() / double(RAND_MAX);
        sys->vx[i] = sin(t1-0.3)/5.0;
        sys->vy[i] = sin(t2-0.3)/5.0;
        sys->vz[i] = sin(t3-0.2)/5.0;
    }
}

void welcome_message()
{
    fprintf(logfile,"Molecular Dynamics Program for simulation fibrinogen polymerization\n");
    fprintf(logfile,"\nUnits used and expected in the data files\n");
    fprintf(logfile,"\n\nDistance -> nm\n");
    fprintf(logfile,"Velocity -> nm/ps\n");
    fprintf(logfile,"Time -> ps\n");
    fprintf(logfile,"Energy -> kJ/mol\n");
    fprintf(logfile,"Mass -> amu or g/mol\n");
    fprintf(logfile,"Charge -> e\n");
    fprintf(logfile,"Force -> kJ/mol/nm\n\n");
    fflush(logfile);
}
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
void init_rcgmd(mdsys_t *sys, rcgmd_t *rcgmd)
{
    int pos, length, type;
    string linebuffer;
    double k1, rbond, k2, rcrit;

    rcgmd->rcrit = 0;
    rcgmd->k1 = 0;
    rcgmd->k2 = 0;
    rcgmd->rbond = 0;
    rcgmd->rcrit_sqr = 0;

    ifstream in;
    in.open("conf.in");
    fprintf(logfile, "\nInitializing rcgmd module..\n");
    //pass 6 gather rcgmd values
    while ( getline(in, linebuffer))
    {
        if(linebuffer.find ("#") != string::npos) continue; // comment line skip it.

        if (linebuffer.find ("pair_sumith") != string::npos)
        {
            //#pair_sumith k1 k2 rbond rcrit
            pos = linebuffer.find(" ");
            linebuffer.erase(0, pos+1); //stripped "pair_sumith"
            //k1 k2 rbond rcrit
            pos = linebuffer.find(" ");
            k1 = atof(linebuffer.substr(0, pos).c_str());
            linebuffer.erase(0, pos+1); //stripped k1
            //k2 rbond rcrit
            pos = linebuffer.find(" ");
            k2 = atof(linebuffer.substr(0, pos).c_str());
            linebuffer.erase(0, pos+1); //stripped k2
            //rbond rcrit
            pos = linebuffer.find(" ");
            rbond = atof(linebuffer.substr(0, pos).c_str());
            linebuffer.erase(0, pos+1); //stripped rbond
            //rcrit
            length = linebuffer.length();
            rcrit = atof(linebuffer.substr(0, length).c_str());
            break;
        }
    }
    in.close();
    rcgmd->rcrit = rcrit;
    rcgmd->k1 = k1;
    rcgmd->k2 = k2;
    rcgmd->rbond = rbond;
    rcgmd->rcrit_sqr =  rcrit * rcrit;
    fprintf(logfile,"rcgmd: k1 = %f, rbond = %f, k2 = %f, d0 = %f\n",rcgmd->k1, rcgmd->rbond, rcgmd->k2, rcgmd->rcrit);

    //create type 2, 3, 4 lists
    rcgmd->ntype2 = 0;
    rcgmd->ntype3 = 0;
    rcgmd->ntype4 = 0;
    for(int i = 0; i < sys->natoms; i++)
    {
        type = sys->type[i];
        if(type == 2) rcgmd->ntype2 += 1;
        if(type == 3) rcgmd->ntype3 += 1;
        if(type == 4) rcgmd->ntype4 += 1;
    }
    int *type2list = (int *) malloc(rcgmd->ntype2 * sizeof(int));
    int *type3list = (int *) malloc(rcgmd->ntype3 * sizeof(int));
    int *type4list = (int *) malloc(rcgmd->ntype4 * sizeof(int));

    int type2 = 0;
    int type3 = 0;
    int type4 = 0;

    for(int i = 0; i < sys->natoms; i++)
    {
        type = sys->type[i];
        if(type == 2)
        {
            type2list[type2] = i;
            type2++;
        }
        if(type == 3)
        {
            type3list[type3] = i;
            type3++;
        }
        if(type == 4)
        {
            type4list[type4] = i;
            type4++;
        }
    }

    rcgmd->type2list = type2list;
    rcgmd->type3list = type3list;
    rcgmd->type4list = type4list;

    //rcgmd->pcgbond = rcgmd->k1 * rcgmd->rbond * rcgmd->rbond;
    //rcgmd->pcg_pair_A = rcgmd->k2 * pow(rcgmd->rcrit, 8.0);
    //rcgmd->pcg_pair_B = rcgmd->k2 * pow(rcgmd->rcrit, 4.0);

    int **neighlist = (int **) malloc(sys->natoms *sizeof(int*));
    for(int i = 0; i < sys->natoms; i++)
        neighlist[i] = (int *) malloc((sys->ntypes+1) * sizeof(int));
    rcgmd->pcgneighlist = neighlist;

    for(int i = 0; i < sys->natoms; i++)
    {
        for(int j = 0; j <=sys->ntypes; j++)
        {
            rcgmd->pcgneighlist[i][j] = 0;
        }
    }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
void rcgmd_force(mdsys_t *sys, rcgmd_t *rcgmd)
{
    int i, j, Nn;
    double rc_sqr = sys->rcut * sys->rcut;
    double rsqr, rx, ry, rz, e_pcg = 0;
    double irr, r, ir, ir4, ir8, ffac;
    double E_new1, E_new2, E_new3, Fnew1, Fnew2, Fnew3;
    double theta1, theta2, theta3, theta4, theta5;
    double theta_block1, theta_block2, theta_block3;
    double rminusrbond, ir_rbon, irr_rbon;

    //-------Phase I----------------------------------
    for(int iatom = 0; iatom < rcgmd->ntype3; iatom++)
    {
        i = rcgmd->type3list[iatom];

        for(int jatom = 0; jatom < rcgmd->ntype4; jatom++)
        {
            j = rcgmd->type4list[jatom];

            rx = pbc(sys->x[i] - sys->x[j], sys->L);
            ry = pbc(sys->y[i] - sys->y[j], sys->B);
            rz = pbc(sys->z[i] - sys->z[j], sys->H);
            rsqr = rx*rx + ry*ry + rz*rz;
            if (rsqr < rc_sqr)
            {
                Nn = max(rcgmd->pcgneighlist[i][4], rcgmd->pcgneighlist[j][3]);
                r = sqrt(rsqr);
                rminusrbond = r - rcgmd->rbond;
                ir_rbon = 1.0 / rminusrbond;
                irr_rbon = ir_rbon * ir_rbon;
                ir = 1.0 / r;
                irr = 1.0 / rsqr;

                E_new1 = -rcgmd->k1 * irr_rbon;
                Fnew1 = -2*rcgmd->k1 * ir * irr_rbon * ir_rbon; // Attractive force

                E_new2 = rcgmd->k2 * rminusrbond * rminusrbond;
                Fnew2 = -2 * rcgmd->k2 * rminusrbond * ir;

                E_new3 = rcgmd->k1 * ir;
                Fnew3 = rcgmd->k1 * ir * irr;  // Repulsive force

                theta1 = heavyside(0.9 - Nn);
                theta2 = heavyside((Nn - 0.9) / (1.1 - Nn));
                theta3 = heavyside(Nn - 1.9);
                theta4 = heavyside(rsqr - rcgmd->rcrit_sqr); // outside region
                theta5 = heavyside(rcgmd->rcrit_sqr - rsqr); // inside region

                //theta_block1 = (theta1 + theta2 * theta5 + theta3 * theta5);
                theta_block1 = theta1 * theta4;
                theta_block2 = theta2 * theta5;
                theta_block3 = (theta3 + theta2 * theta4);

                //fprintf(logfile, "%d\t%d\t%d,%d\t%d\t%d\t%d\t%d\n",theta_block1, theta_block2,theta_block3, theta1,theta2,theta3,theta4,theta5);

                e_pcg += theta_block1 * E_new1 + theta_block2 * E_new2 + theta_block3 * E_new3;
                ffac = theta_block1 * Fnew1 + theta_block2 * Fnew2 + theta_block3 * Fnew3;

                sys->fx[i] += rx * ffac;
                sys->fy[i] += ry * ffac;
                sys->fz[i] += rz * ffac;

                sys->fx[j] -= rx * ffac;
                sys->fy[j] -= ry * ffac;
                sys->fz[j] -= rz * ffac;
            }
        }
    }
    sys->epot += e_pcg;

    //-------Phase II----------------------------------
    e_pcg = 0;
    for(int iatom = 0; iatom < rcgmd->ntype2; iatom++)
    {
        i = rcgmd->type2list[iatom];

        for(int jatom = 0; jatom < rcgmd->ntype4; jatom++)
        {
            j = rcgmd->type4list[jatom];

            if(rcgmd->pcgneighlist[j][3] == 1) //# work only thrombins (4) attached to E region (3)
            {
                int Nn24, Nn42;

                rx = pbc(sys->x[i] - sys->x[j], sys->L);
                ry = pbc(sys->y[i] - sys->y[j], sys->B);
                rz = pbc(sys->z[i] - sys->z[j], sys->H);
                rsqr = rx*rx + ry*ry + rz*rz;
                if (rsqr < rc_sqr)
                {
                    Nn24 = rcgmd->pcgneighlist[i][4];  //# type 4 neighbors around i (2)
                    Nn42 = rcgmd->pcgneighlist[j][2];  //# type 2 neighbors around j (4)
                    // Nn42 is always 1 here---optimize thetas
                    r = sqrt(rsqr);
                    rminusrbond = r - rcgmd->rbond;
                    ir_rbon = 1.0 / rminusrbond;
                    irr_rbon = ir_rbon * ir_rbon;
                    ir = 1.0 / r;
                    irr = 1.0 / rsqr;

                    E_new1 = -rcgmd->k1 * irr_rbon;
                    Fnew1 = -2*rcgmd->k1 * ir * irr_rbon * ir_rbon; // Attractive force

                    E_new2 = rcgmd->k2 * rminusrbond * rminusrbond;
                    Fnew2 = -2 * rcgmd->k2 * rminusrbond * ir;

                    E_new3 = rcgmd->k1 * ir;
                    Fnew3 = rcgmd->k1 * ir * irr;  // Repulsive force

                    theta1 = heavyside(0.9 - Nn42) * heavyside(0.9 - Nn24);
                    theta2 = (heavyside((Nn42 - 0.9) / (1.1 - Nn42)) + heavyside((Nn42 - 1.9) / (2.1 - Nn42))) * heavyside((Nn24 - 0.9) / (1.1 - Nn24));
                    theta3 = max(heavyside(Nn24 - 1.9), heavyside(Nn42 - 2.9));
                    theta4 = heavyside(rsqr - rcgmd->rcrit_sqr);
                    theta5 = heavyside(rcgmd->rcrit_sqr - rsqr);

                    //theta_block1 = (theta1 + theta2 * theta5 + theta3 * theta5);
                    theta_block1 = theta1 * theta4;
                    theta_block2 = theta2 * theta5;
                    theta_block3 = (theta3 + theta2 * theta4);

                    //fprintf(logfile, "%d\t%d\t%d,%d\t%d\t%d\t%d\t%d\n",theta_block1, theta_block2,theta_block3, theta1,theta2,theta3,theta4,theta5);

                    e_pcg += theta_block1 * E_new1 + theta_block2 * E_new2 + theta_block3 * E_new3;
                    ffac = theta_block1 * Fnew1 + theta_block2 * Fnew2 + theta_block3 * Fnew3;

                    sys->fx[i] += rx * ffac;
                    sys->fy[i] += ry * ffac;
                    sys->fz[i] += rz * ffac;

                    sys->fx[j] -= rx * ffac;
                    sys->fy[j] -= ry * ffac;
                    sys->fz[j] -= rz * ffac;
                }
            }
        }
    }
    sys->epot += e_pcg;
    //-------Phase III----------------------------------

}
//----------------------------------------------------------------------
double heavyside(double val)
{
    if(val < 0)    return 0;
    else if(val == 0)  return 0.5;
    else    return 1;
}
