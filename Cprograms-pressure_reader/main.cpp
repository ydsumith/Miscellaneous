#include <conio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <string.h>
#include <stdio.h>
#include <cmath>

using namespace std;


int error_status = 0;
int NATOMS, NZ;
int N_FRAMES = 0;

double xlo, xhi, ylo, yhi, zlo, zhi;
double * PXX;
double * PYY;
double * PZZ;
double * PXY;
double * PYZ;
double * PXZ;
double * NDEN;
double dz = 1; // 1 Ang is the slab thickness
double factor = 100000; //divide by this number to avoid possible large numbers


string FILENAME = "dump_press.mol";
ifstream in;
ofstream output_file;

int file_read_pass1();
int construct_memory();
int destruct_memory();
int read_and_summation();
int print_results();
//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------
int main()
{
    cout << "----------Welcome to the pressure estimator from MD simulations--------" << endl;
    cout << "[IMPORTANT] Please make sure you have dump_press.mol in the current directory of this exe file"<<endl;
    cout << "\nAre you sure to continue? (y/n): "<<endl;
    if(getch()== 'y')
    {
        cout<<"Please be patient while I get through this mess..."<<endl;
        error_status = file_read_pass1();
        if (error_status == 1)
        {
            cout<<"[ERROR] File read error."<<endl;
        }
        construct_memory();
        read_and_summation();
        print_results();
        cout<<"Program completed."<<endl;
    }
    else
    {
        cout<<"Exiting as per your wish.."<<endl;
    }

    destruct_memory();
    return 0;
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
int print_results()
{
    output_file.open ("results_unavg.txt");
    output_file << "dz\tPXX\tPYY\tPZZ\tPXY\tPXZ\tPYZ\tNDEN\n";
    for(int i = 0; i < NZ ; i++)
    {
        if(NDEN[i]!= 0)
        {
            PXX[i] = PXX[i] * factor / NDEN[i];
            PYY[i] = PYY[i] * factor / NDEN[i];
            PZZ[i] = PZZ[i] * factor / NDEN[i];
            PXY[i] = PXY[i] * factor / NDEN[i];
            PXZ[i] = PXZ[i] * factor / NDEN[i];
            PYZ[i] = PYZ[i] * factor / NDEN[i];
        }
    }

    for(int i = 0; i < NZ ; i++)
    {
        output_file<<i*dz<<"\t"<<PXX[i]<<"\t"<<PYY[i]<<"\t"<<PZZ[i]<<"\t"<<PXY[i]<<"\t"<<PXZ[i]<<"\t"<<PYZ[i]<<"\t"<<NDEN[i]<<"\n";
    }
    output_file.close();

	//---------------------------------------------------------
	output_file.open ("results_avg.txt");
    output_file << "dz\tPXX\tPYY\tPZZ\tPXY\tPXZ\tPYZ\n";
    for(int i = 0; i < NZ ; i++)
    {
		PXX[i] = PXX[i] / N_FRAMES;
		PYY[i] = PYY[i] / N_FRAMES;
		PZZ[i] = PZZ[i] / N_FRAMES;
		PXY[i] = PXY[i] / N_FRAMES;
		PXZ[i] = PXZ[i] / N_FRAMES;
		PYZ[i] = PYZ[i] / N_FRAMES;
    }

    for(int i = 0; i < NZ ; i++)
    {
		output_file<<i*dz<<"\t"<<PXX[i]<<"\t"<<PYY[i]<<"\t"<<PZZ[i]<<"\t"<<PXY[i]<<"\t"<<PXZ[i]<<"\t"<<PYZ[i]<<"\n";
    }
    output_file.close();
    return 0;
}

//--------------------------------------------------------------------------
int read_and_summation()
{
    int counter = 0;
    int arry_cntr;
    int length, startLOC, endLOC, zLOC;
    int BLK_LINES = 9+NATOMS;
    double temp_holder[15];
    string linebuffer;

    in.open(FILENAME.c_str());
    while ( getline(in, linebuffer))
    {
        counter++;
        if((counter % BLK_LINES) > 9) // pressure information starts here
        {
            length = linebuffer.length();
            startLOC = endLOC = 0;
            for(int i = 0, arry_cntr = 0; i < length; i++)
            {
                if(linebuffer.c_str()[i]== ' ')
                {
                    endLOC = i;
                    temp_holder[arry_cntr] = atof(linebuffer.substr(startLOC, endLOC).c_str());
                    arry_cntr++;
                    startLOC = i+1;
                }
            }
            zLOC = floor((temp_holder[5]-zlo)/(dz)); // getting the z location
            //xx, yy, zz, xy, xz, yz.
            PXX[zLOC] += temp_holder[6]/factor;
            PYY[zLOC] += temp_holder[7]/factor;
            PZZ[zLOC] += temp_holder[8]/factor;
            PXY[zLOC] += temp_holder[9]/factor;
            PXZ[zLOC] += temp_holder[10]/factor;
            PYZ[zLOC] += temp_holder[11]/factor;
            NDEN[zLOC] += 1;
        }

        if((counter % BLK_LINES) == 0)
		{
			N_FRAMES += 1;
			//cout<<"N_FRAMES = "<<N_FRAMES<<endl;
		}
    }
    in.close();
    return 0;
}

//--------------------------------------------------------------------------
int construct_memory()
{
    PXX = new (nothrow) double [NZ];
    PYY = new (nothrow) double [NZ];
    PZZ = new (nothrow) double [NZ];
    PXY = new (nothrow) double [NZ];
    PYZ = new (nothrow) double [NZ];
    PXZ = new (nothrow) double [NZ];
    NDEN = new (nothrow) double [NZ];

    for(int i = 0; i < NZ; i++)
    {
        PXX[i] = 0;
        PYY[i] = 0;
        PZZ[i] = 0;
        PXY[i] = 0;
        PYZ[i] = 0;
        PXZ[i] = 0;
        NDEN[i] = 0;
    }
    return 0;
}

//--------------------------------------------------------------------------
int destruct_memory()
{
    delete PXX;
    delete PYY;
    delete PZZ;
    delete PXY;
    delete PYZ;
    delete PXZ;
    delete NDEN;
    return 0;
}

//--------------------------------------------------------------------------
int file_read_pass1()
{
    int counter = 0, length, pos;
    string linebuffer, delim;

    in.open(FILENAME.c_str());

    while ( getline(in, linebuffer))
    {
        counter++;
        if(counter > 8)
            break;
        if(counter == 4) //no: of atoms
        {
            length = linebuffer.length();
            NATOMS = atoi(linebuffer.substr(0, length).c_str());
            cout<<counter<<","<<NATOMS<<","<<linebuffer<<endl;
        }
        if(counter == 6) // xlo xhi
        {
            length = linebuffer.length();
            if (linebuffer.find (" ") != string::npos)
                delim = " ";
            else if (linebuffer.find ("\t") != string::npos)
                delim = "\t";

            pos = linebuffer.find(delim);

            xlo = atof(linebuffer.substr(0, pos).c_str());
            xhi = atof(linebuffer.substr(pos+1, length).c_str());
            cout<<counter<<","<<xlo<<","<<xhi<<endl;
        }
        if(counter == 7) // ylo yhi
        {
            length = linebuffer.length();
            if (linebuffer.find (" ") != string::npos)
                delim = " ";
            else if (linebuffer.find ("\t") != string::npos)
                delim = "\t";

            pos = linebuffer.find(delim);

            ylo = atof(linebuffer.substr(0, pos).c_str());
            yhi = atof(linebuffer.substr(pos+1, length).c_str());
            cout<<counter<<","<<ylo<<","<<yhi<<endl;
        }
        if(counter == 8) // zlo zhi
        {
            length = linebuffer.length();
            if (linebuffer.find (" ") != string::npos)
                delim = " ";
            else if (linebuffer.find ("\t") != string::npos)
                delim = "\t";

            pos = linebuffer.find(delim);

            zlo = atof(linebuffer.substr(0, pos).c_str());
            zhi = atof(linebuffer.substr(pos+1, length).c_str());
            cout<<counter<<","<<zlo<<","<<zhi<<endl;
        }
    }
    in.close();

    NZ = int((zhi-zlo)/dz);
    cout<<"NZ = "<<NZ<<endl;
    return 0;
}
