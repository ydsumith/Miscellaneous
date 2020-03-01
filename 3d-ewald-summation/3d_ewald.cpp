/********************************************************
 * All formulae taken from 
 * Essmann, Ulrich, et al. "A smooth particle mesh Ewald method." 
 * The Journal of chemical physics 103.19 (1995): 8577-8593.
********************************************************/
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <cstring>
#include <string.h>
#include <iomanip>
#include <omp.h>


using namespace std;

ofstream output;

int VAL, N;

double U_CR, U_CF, U_CS, U_CLR, U_INTRA;
double Vir_XX, Vir_YY, Vir_ZZ;
double PXX, PYY, PZZ;
double L,B,H,ALPHA;
double ONE_PI_EPS0, half_L, half_B, half_H;
double nmax, kmax, box_VOL;

bool build_geometry();
void define_physics();
void REAL_SPACE();
void FOURIER_SPACE();
void SELF();
void DIPOLE();
void DISP_FORCE();
void read_water();
void INTRA_MOL();
void get_PBC(double &RX, double &RY, double &RZ);
void find_COM();
void virial_com();

struct orion
{
    double x;
    double y;
    double z;
    double charge;
    double mass;
    double C12_root;
    double C6_root;
    int molecule;
};
struct nebula
{
    double x;
    double y;
    double z;
};

//----------declaration of independence--------
nebula *R_force, *K_force, *Intra_force, *FORCE, *DIP_force;
nebula *COM;
orion *atoms;
//------------------------------------------------------------------
//
//            The legendary Main function
//
//------------------------------------------------------------------
int main()
{
    output<<"/******* Ewald summation program******/\n Authored by Sumith YD\n";

    define_physics();
    // Get the number of processors in this system
    int iCPU = omp_get_num_procs();

    // Now set the number of threads
    omp_set_num_threads(iCPU); //--- not relevant now. may be for future.

    if(build_geometry()) output<<"Geometry has builded"<<endl; // read gro file with coordinates

    find_COM(); // find center of mass (optional)

    REAL_SPACE(); //real space contributions
    output<<"REAL ENERGY = "<<U_CR<<" [kJ/mol]"<<endl;

    FOURIER_SPACE(); //reciprocal space
    output<<"K_SPACE ENERGY = "<<U_CF<<" [kJ/mol]"<<endl;

    SELF(); //self energy
    output<<"SELF ENERGY = "<<U_CS<<" [kJ/mol]"<<endl;

    DIPOLE(); //dipole contributions (if any)
    output<<"DIPOLE ENERGY = "<<U_CLR<<" [kJ/mol]"<<endl<<endl;

    INTRA_MOL(); // correction for intra molecular forces
    output<<"INTRA ENERGY = "<<U_INTRA<<" [kJ/mol]"<<endl<<endl;
    output<<"K_SPACE+SELF = "<<U_CF+U_CS<<" [kJ/mol]"<<endl;
    output<<"TOTAL COUL ENERGY = "<<U_CF+U_CS+U_CR+U_INTRA<<" [kJ/mol]"<<endl;

    DISP_FORCE(); // write results to a file
    output<<"Mission accomplished.."<<endl<<endl<<endl;
    output.close();
    return 0;
}

//---------------------------------------------------------------------------
//
//
//---------------------------------------------------------------------------
void REAL_SPACE()
{
    double qi, qj, rsqr, r;
    double rx,ry,rz, tmp;
    double RX, RY, RZ;
    double GAMMA = 2*ALPHA/sqrt(M_PI);
    double ALP2 = ALPHA*ALPHA;

    U_CR = 0;
    output<<"calculating real space values..."<<endl;
    for(int i = 0; i < N-1; i++)
    {
        qi = atoms[i].charge;
        for(int j = i+1; j < N; j++)
        {
            qj = atoms[j].charge;
            rx = atoms[i].x - atoms[j].x;
            ry = atoms[i].y - atoms[j].y;
            rz = atoms[i].z - atoms[j].z;
            get_PBC(rx, ry, rz);
            for(int nx = -nmax; nx <= nmax ; nx++)
            {
                for(int ny = -nmax; ny <= nmax ; ny++)
                {
                    for(int nz = -nmax; nz <= nmax ; nz++)
                    {
                        if(nx == 0 && ny == 0 && nz == 0) //central cell
                        {
                            //ignore i == j
                            if(atoms[i].molecule != atoms[j].molecule)
                            {
                                rsqr = rx*rx + ry*ry + rz*rz;
                                if(rsqr != 0)
                                {
                                    r = sqrt(rsqr);
                                    U_CR += ONE_PI_EPS0 * qi * qj * erfc(ALPHA * r) / r;
                                    tmp = (GAMMA * exp(-ALP2*rsqr) + erfc(ALPHA*r)/r)/rsqr;
                                    tmp = tmp * ONE_PI_EPS0 * qi * qj;
                                    R_force[i].x += tmp * rx;
                                    R_force[i].y += tmp * ry;
                                    R_force[i].z += tmp * rz;

                                    R_force[j].x -= tmp * rx;
                                    R_force[j].y -= tmp * ry;
                                    R_force[j].z -= tmp * rz;

                                    /********* virial *********/
                                    Vir_XX += tmp * rx* rx;
                                    Vir_YY += tmp * ry* ry;
                                    Vir_ZZ += tmp * rz* rz;
                                }
                            }
                            else
                            {
                                //ignore, because same molecule
                            }
                        }
                        else //image cells
                        {
                            RX = rx + nx*L;
                            RY = ry + ny*B;
                            RZ = rz + nz*H;
                            rsqr = RX*RX + RY*RY + RZ*RZ;
                            if(rsqr != 0)
                            {
                                r = sqrt(rsqr);
                                U_CR += ONE_PI_EPS0 * qi * qj * erfc(ALPHA * r) / r;
                                tmp = (GAMMA * exp(-ALP2*rsqr) + erfc(ALPHA*r)/r)/rsqr;
                                tmp = tmp * ONE_PI_EPS0 * qi * qj;
                                R_force[i].x += tmp * RX;
                                R_force[i].y += tmp * RY;
                                R_force[i].z += tmp * RZ;

                                R_force[j].x -= tmp * RX;
                                R_force[j].y -= tmp * RY;
                                R_force[j].z -= tmp * RZ;

                                /********* virial *********/
                                Vir_XX += tmp * RX* RX;
                                Vir_YY += tmp * RY* RY;
                                Vir_ZZ += tmp * RZ* RZ;
                            }
                        }
                    }
                }
            }
        }
    }
    output<<"done with real space"<<endl;
}
//---------------------------------------------------------------------------
//
//
//---------------------------------------------------------------------------
void FOURIER_SPACE()
{
    complex<double> ak;
    double ksqr, kdotr, qi, tmp, tmp2;
    double a,b,akak, mx,my,mz, chak;
    double GAMMA = -1/(4*ALPHA*ALPHA);
    double recip = ONE_PI_EPS0 * 2* M_PI / box_VOL;
    //output.open("DEBUG_details.xls");
    //output<<"i\tkz\tmz\tqi\tkdotr\ttmp2\tforcez"<<endl;
    U_CF = 0;

    for(int kx = -kmax; kx <= kmax ; kx++)
    {
        mx = 2*M_PI*kx/L;
        for(int ky = -kmax; ky <= kmax ; ky++)
        {
            my = 2*M_PI*ky/B;
            for(int kz = -kmax; kz <= kmax ; kz++)
            {
                mz = 2*M_PI*kz/H;
                ksqr = mx*mx + my*my + mz*mz;
                if(ksqr != 0)
                {
                    ak.real() = 0;
                    ak.imag() = 0;
                    for(int i = 0; i < N; i++)
                    {
                        kdotr = mx*atoms[i].x + my*atoms[i].y + mz*atoms[i].z;
                        qi = atoms[i].charge;
                        ak.real() += qi*cos(kdotr);
                        ak.imag() -= qi*sin(kdotr);
                    }
                    a = ak.real();
                    b = ak.imag();
                    akak = (a*a + b*b);
                    tmp = recip * exp(GAMMA * ksqr)/ksqr;
                    U_CF += tmp * akak;
                    chak = (2/ksqr) - 2*GAMMA;
                    /********* virial *********/
                    Vir_XX += (1 - chak * mx*mx) *tmp * akak;
                    Vir_YY += (1 - chak * my*my) *tmp * akak;
                    Vir_ZZ += (1 - chak * mz*mz) *tmp * akak;
                    for(int i = 0; i < N; i++)
                    {
                        kdotr = mx*atoms[i].x + my*atoms[i].y + mz*atoms[i].z;
                        qi = atoms[i].charge;
                        tmp2 = 2*tmp*qi*(sin(kdotr) * a + cos(kdotr) * b);
                        K_force[i].x += tmp2 * mx;
                        K_force[i].y += tmp2 * my;
                        K_force[i].z += tmp2 * mz;
                        //output<<setprecision(15)<<i<<"\t"<<kz<<"\t"<<mz<<"\t"<<qi<<"\t"<<kdotr<<"\t"<<tmp2<<"\t"<<tmp2*mz<<"\t"<<K_force[i].z<<endl;
                    }
                }
            }
        }
    }
    //output.close();
}
//---------------------------------------------------------------------------
//
//
//---------------------------------------------------------------------------
void DIPOLE()
{
    double M_sqr, Mx, My, Mz;
    double qi ;

    U_CLR = 0;
    M_sqr = Mx = My = Mz = 0;
    for(int i = 0; i < N; i++)
    {
        qi = atoms[i].charge;
        Mx += qi * atoms[i].x;
        My += qi * atoms[i].y;
        Mz += qi * atoms[i].z;
    }
    M_sqr = Mx*Mx + My*My + Mz*Mz;
    U_CLR = ONE_PI_EPS0 * 2 * M_PI * M_sqr / (box_VOL*3);
    for(int i = 0; i < N; i++)
    {
        qi = atoms[i].charge;
        DIP_force[i].x -= ONE_PI_EPS0* qi * 4*M_PI * Mx / (3*box_VOL);
        DIP_force[i].y -= ONE_PI_EPS0* qi * 4*M_PI * My / (3*box_VOL);
        DIP_force[i].z -= ONE_PI_EPS0* qi * 4*M_PI * Mz / (3*box_VOL);
        /********* virial *********/
        Vir_XX += ONE_PI_EPS0 * 2*M_PI*(M_sqr - 2*Mx*Mx)/(box_VOL*box_VOL*3);
        Vir_YY += ONE_PI_EPS0 * 2*M_PI*(M_sqr - 2*My*My)/(box_VOL*box_VOL*3);
        Vir_ZZ += ONE_PI_EPS0 * 2*M_PI*(M_sqr - 2*Mz*Mz)/(box_VOL*box_VOL*3);
    }
}
//---------------------------------------------------------------------------
//
//
//---------------------------------------------------------------------------
void SELF()
{
    double SR_PI = -ALPHA/sqrt(M_PI);

    U_CS = 0;
    for(int i = 0; i < N; i++)
    {
        U_CS += atoms[i].charge * atoms[i].charge;
    }
    U_CS *= SR_PI * ONE_PI_EPS0;
}
//---------------------------------------------------------------------------
//
//
//---------------------------------------------------------------------------
void INTRA_MOL()
{
    double rx, ry, rz, qi, qj;
    double rsqr, r, tmp;
    double KAPPA = 2*ALPHA/sqrt(M_PI);
    double ALPSQR = ALPHA * ALPHA;
    double FST;

    U_INTRA = 0;

    for(int i = 0; i < N-1; i++) // N is the number of atoms
    {
        for(int j = i+1; j < N; j++)
        {
            if(atoms[i].molecule == atoms[j].molecule) // same molecule
            {
                qi = atoms[i].charge;
                qj = atoms[j].charge;
                rx = atoms[i].x - atoms[j].x;
                ry = atoms[i].y - atoms[j].y;
                rz = atoms[i].z - atoms[j].z;
                rsqr = rx*rx + ry*ry + rz*rz;
                r = sqrt(rsqr);
                FST =erf(ALPHA*r);
                U_INTRA -= ONE_PI_EPS0*qi*qj*FST/r;
                tmp = ONE_PI_EPS0*qi*qj*(KAPPA*exp(-ALPSQR*rsqr)/rsqr - FST/(r*rsqr));

                Intra_force[i].x += tmp * rx;
                Intra_force[i].y += tmp * ry;
                Intra_force[i].z += tmp * rz;

                Intra_force[j].x -= tmp * rx;
                Intra_force[j].y -= tmp * ry;
                Intra_force[j].z -= tmp * rz;
                /********* virial *********/
                Vir_XX += tmp * rx * rx;
                Vir_YY += tmp * ry * ry;
                Vir_ZZ += tmp * rz * rz;
            }
        }
    }
}
//---------------------------------------------------------------------------
//
//
//---------------------------------------------------------------------------
void define_physics()
{
    ONE_PI_EPS0 =138.9354859;
}
//---------------------------------------------------------------------------
//
//
//---------------------------------------------------------------------------
bool build_geometry()
{
    output.open("LOGFILE.txt");
    nmax = 5;
    kmax = 13; //update for future

    read_water();

    ALPHA = 1/0.313759; // this has to be optimized based on cutoff
    box_VOL = L*B*H;

    for(int i = 0; i < N; i++)
    {
        R_force[i].x = 0;
        R_force[i].y = 0;
        R_force[i].z = 0;
        K_force[i].x = 0;
        K_force[i].y = 0;
        K_force[i].z = 0;
        Intra_force[i].x = 0;
        Intra_force[i].y = 0;
        Intra_force[i].z = 0;
        DIP_force[i].x = 0;
        DIP_force[i].y = 0;
        DIP_force[i].z = 0;
        FORCE[i].x = 0;
        FORCE[i].y = 0;
        FORCE[i].z = 0;
    }
    PXX = PYY = PZZ = 0;
    Vir_XX = Vir_YY = Vir_ZZ = 0;
    return true;
}
//---------------------------------------------------------------------------
//
//
//---------------------------------------------------------------------------
void read_water()
{
    string str_val,atom_type;
    int i = 0, fnd;
    size_t found;
    ifstream readit;

    readit.open("conf.gro",ios::in);
    getline(readit,str_val); //reading first line and ignore
    getline(readit,str_val); //reading the total # of atoms
    N = atoi(str_val.c_str());
    R_force = new (nothrow) nebula[N];
    K_force = new (nothrow) nebula[N];
    Intra_force = new (nothrow) nebula[N];
    DIP_force = new (nothrow) nebula[N];
    FORCE = new (nothrow) nebula[N];
    atoms = new (nothrow) orion[N];
    COM = new (nothrow) nebula[N];

    while(getline(readit,str_val))
    {
        if(i < N)
        {
            atom_type = str_val.substr(11,5).c_str();
            fnd = 0;
            found = atom_type.find("OW");
            if(found != string::npos)
            {
                atoms[i].charge = -0.8476;
                atoms[i].C12_root = sqrt(4*0.65*pow(0.3166,12));
                atoms[i].C6_root = sqrt(4*0.65*pow(0.3166,6));
                atoms[i].mass = 15.99491461956;
                fnd =1;
            }
            found = atom_type.find("HW");
            if(found != string::npos)
            {
                atoms[i].charge = 0.4238;
                atoms[i].C12_root = atoms[i].C6_root = 0;
                atoms[i].mass = 1.0078250;
                fnd =1;
            }
            if(fnd == 0)
            {
                output<<"[ERROR] Corrupted gro file..."<<endl<<"exiting...";
                output<<"N = "<<N<<"\n\n";
                exit(0);
            }
            atoms[i].molecule = atoi(str_val.substr(1,5).c_str());
            atoms[i].x = atof(str_val.substr(21,8).c_str());
            atoms[i].y = atof(str_val.substr(29,8).c_str());
            atoms[i].z = atof(str_val.substr(37,8).c_str());
            //md_log<<i<<","<<atoms[i].x  <<","<< atoms[i].y <<","<<atoms[i].z <<endl;
            i++;
        }
        //code for L,B,H here
        L = atof(str_val.substr(1,10).c_str());
        B = atof(str_val.substr(11,10).c_str());
        H = atof(str_val.substr(21,10).c_str());
    }
    readit.close();
    half_L = L/2;
    half_B = B/2;
    half_H = H/2;
}
//---------------------------------------------------------------------------
//
//
//---------------------------------------------------------------------------
void get_PBC(double &RX, double &RY, double &RZ)
{
    if (abs(RX) > half_L)
    {
        if(RX <0 )	RX = RX + L;
        else		RX = RX - L;
    }
    if (abs(RY) > half_B)
    {
        if(RY <0 )	RY = RY + B;
        else		RY = RY - B;
    }
    if (abs(RZ) > half_H)
    {
        if(RZ <0 )	RZ = RZ + H;
        else		RZ = RZ - H;
    }
}
//---------------------------------------------------------------------------
//
//
//---------------------------------------------------------------------------
void find_COM()
{
    double m1, m2, m3;
    for(int i = 0; i < N; i += 3)
    {
        m1 = atoms[i].mass;
        m2 = atoms[i+1].mass;
        m3 = atoms[i+2].mass;
        COM[i].x = (atoms[i].x*m1 +atoms[i+1].x*m2 +atoms[i+2].x*m3)/(m1+m2+m3);
        COM[i].y = (atoms[i].y*m1 +atoms[i+1].y*m2 +atoms[i+2].y*m3)/(m1+m2+m3);
        COM[i].z = (atoms[i].z*m1 +atoms[i+1].z*m2 +atoms[i+2].z*m3)/(m1+m2+m3);
        COM[i+1].x = COM[i+2].x = COM[i].x;
        COM[i+1].y = COM[i+2].y = COM[i].y;
        COM[i+1].z = COM[i+2].z = COM[i].z;
    }
}
//---------------------------------------------------------------------------
//
//
//---------------------------------------------------------------------------
void virial_com() //virial correction for rigid molecules
{
    for(int i = 0; i < N; i++)
    {
        Vir_XX -= FORCE[i].x *(atoms[i].x - COM[i].x);
        Vir_YY -= FORCE[i].y *(atoms[i].y - COM[i].y);
        Vir_ZZ -= FORCE[i].z *(atoms[i].z - COM[i].z);
    }
}
//---------------------------------------------------------------------------
//
//
//---------------------------------------------------------------------------
void DISP_FORCE()
{
    output<<endl<<"FORCE [kJ/mol/nm]"<<endl;
    output<<"atom\tFX\tFY\tFZ "<<endl;
    for(int i = 0; i < N; i++)
    {
        FORCE[i].x = R_force[i].x + K_force[i].x + Intra_force[i].x;
        FORCE[i].y = R_force[i].y + K_force[i].y + Intra_force[i].y;
        FORCE[i].z = R_force[i].z + K_force[i].z + Intra_force[i].z;
    }
    for(int i = 0; i < N; i++)
        output<<i+1<<"\t"<<FORCE[i].x<<"\t"<<FORCE[i].y<<"\t"<<FORCE[i].z<<endl;

    output<<endl<<"R_FORCE [kJ/mol/nm]"<<endl;
    output<<"atom\tFX\tFY\tFZ "<<endl;
    for(int i = 0; i < N; i++)
        output<<i+1<<"\t"<<R_force[i].x<<"\t"<<R_force[i].y<<"\t"<<R_force[i].z<<endl;

    output<<endl<<"K_FORCE [kJ/mol/nm]"<<endl;
    output<<"atom\tFX\tFY\tFZ "<<endl;

    for(int i = 0; i < N; i++)
        output<<i+1<<"\t"<<K_force[i].x<<"\t"<<K_force[i].y<<"\t"<<K_force[i].z<<endl;

    output<<endl<<"INT_FORCE [kJ/mol/nm]"<<endl;
    output<<"atom\tFX\tFY\tFZ "<<endl;
    for(int i = 0; i < N; i++)
        output<<i+1<<"\t"<<Intra_force[i].x<<"\t"<<Intra_force[i].y<<"\t"<<Intra_force[i].z<<endl;

    output<<endl<<"DIPOL_FORCE [kJ/mol/nm]"<<endl;
    output<<"atom\tFX\tFY\tFZ "<<endl;
    for(int i = 0; i < N; i++)
        output<<i+1<<"\t"<<DIP_force[i].x<<"\t"<<DIP_force[i].y<<"\t"<<DIP_force[i].z<<endl;
    /********* virial *********/
    virial_com();
    output<<endl;
    output<<endl;
    output<<endl<<"VIRIAL [kJ/mol]"<<endl;
    output<<"Vir-XX\tVir-YY\tVir-ZZ"<<endl;
    output<<Vir_XX<<"\t"<<Vir_YY<<"\t"<<Vir_ZZ<<endl;
    PXX = 16.6054 * Vir_XX / box_VOL;
    PYY = 16.6054 * Vir_YY / box_VOL;
    PZZ = 16.6054 * Vir_ZZ / box_VOL;
    output<<endl<<"PRESSURE [bar]"<<endl;
    output<<"PXX\tPYY\tPZZ"<<endl;
    output<<PXX<<"\t"<<PYY<<"\t"<<PZZ<<endl;
}
