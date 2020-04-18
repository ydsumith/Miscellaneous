"""
Created on Fri Apr 17 18:40:38 2020
@author: Prof. Sumith Yesudasan
"""
import numpy as np
import matplotlib.pyplot as plt
import math


def main():
    print("\ntest with values close to xlo, xhi, zlo, zhi")
    print("check sum == 1\n")
    xlo = -5
    xhi = 5
    zlo = -3
    zhi = 3
    dr = 0.5

    NX = int((xhi-xlo)/dr)+2
    NZ = int((zhi-zlo)/dr)+2

    x = np.linspace(xlo, xhi, NX)
    y = np.linspace(zlo, zhi, NZ)

    GX, GY = np.meshgrid(x, y)
    GZ = np.zeros([NZ, NX])

    datax = -4.985
    datay = -2.999

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

    print(sum(sum(GZ)))

    plt.contourf(GX, GY, GZ, 20, cmap='binary')
    plt.colorbar()

    print("Finished the program")


if __name__ == "__main__":
    main()
