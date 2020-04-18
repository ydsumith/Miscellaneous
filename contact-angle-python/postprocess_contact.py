import os
import numpy as np
import re
import matplotlib.pyplot as plt
from scipy.interpolate import Rbf
import math

def get_current_folder():
    cur_dir1 = os.getcwd()  # works with Spyder
    cur_dir2 = cur_dir1.replace("\\", "/")
    return cur_dir1, cur_dir2


def GetColour(v,  vmin,  vmax):
    r = 1
    g = 1
    b = 1
    if (v < vmin):
        v = vmin
    if (v > vmax):
        v = vmax

    dv = vmax - vmin
    if (v < (vmin + 0.25 * dv)):
        r = 0
        g = 4 * (v - vmin) / dv
    elif (v < (vmin + 0.5 * dv)):
        r = 0
        b = 1 + 4 * (vmin + 0.25 * dv - v) / dv
    elif (v < (vmin + 0.75 * dv)):
        r = 4 * (v - vmin - 0.5 * dv) / dv
        b = 0
    else:
        g = 1 + 4 * (vmin + 0.75 * dv - v) / dv
        b = 0
    return r, g, b


def euclidean_norm_numpy(x1, x2):
    return np.linalg.norm(x1 - x2, axis=0)


def main():
    epsilon = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
                        0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    zlo = 0
    zhi = 200
    xlo = -75
    xhi = 75
    dr = 0.5

    NX = int((xhi-xlo)/dr)+5
    NZ = int((zhi-zlo)/dr)+5    
    x = np.linspace(xlo, xhi, NX)
    y = np.linspace(zlo, zhi, NZ)
    GX, GY = np.meshgrid(x, y)

    dir1, dir2 = get_current_folder()

    for Temp in epsilon:
        in_file = dir2 + "/data/dump_" + str(Temp) + ".xyz"
        out_file = dir2 + "/results/result_" + str(Temp) + ".txt"
        out_image = dir2 + "/results/img_interp_" + str(Temp) + ".png"
        print("File to read is \n", in_file)
        print("File to write is \n", out_file)
        
        GZ = np.zeros([NZ, NX])
        
        print("\nPlease wait while I am working...")
        with open(in_file) as file_in:
            for line in file_in:
                temp1 = re.findall(
                    r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)
                # An expected line should be 2 -47.25 -74.25 0 [type x y z]
                if len(temp1) == 4:  # for a valid entry
                    if float(temp1[0]) == 1:  # for type 1, water beads
                        datax = float(temp1[1])
                        datay = float(temp1[3])
                        xg2 = math.floor((datax-xlo)/dr)+1
                        yg2 = math.floor((datay-zlo)/dr)+1
                        xg1 = xg2 - 1
                        yg1 = yg2 - 1
                        xg3 = xg2 + 1
                        yg3 = yg2 + 1
                    
                        xxp = xg2*dr+xlo - datax
                        w1x = 0.125*(1-4*xxp+4*xxp*xxp)
                        w2x = 0.25*(3-4*xxp*xxp)
                        w3x = 0.125*(1+4*xxp+4*xxp*xxp)
                        
                        yyp = yg2*dr+zlo - datay
                        w1y = 0.125*(1-4*yyp+4*yyp*yyp)
                        w2y = 0.25*(3-4*yyp*yyp)
                        w3y = 0.125*(1+4*yyp+4*yyp*yyp)
                        
                        GZ[xg1][yg1] += w1x*w1y
                        GZ[xg1][yg2] += w1x*w2y
                        GZ[xg1][yg3] += w1x*w3y
                        
                        GZ[xg2][yg1] += w2x*w1y
                        GZ[xg2][yg2] += w2x*w2y
                        GZ[xg2][yg3] += w2x*w3y
                        
                        GZ[xg3][yg1] += w3x*w1y
                        GZ[xg3][yg2] += w3x*w2y
                        GZ[xg3][yg3] += w3x*w3y

        # ---plotting business
        plt.show()
        plt.contourf(GX, GY, GZ, 20, cmap='binary')
        plt.colorbar()
        plt.axis(aspect='image')
        plt.savefig(out_image, dpi=1200)

    print("\nProgram finished.")


if __name__ == "__main__":
    main()
