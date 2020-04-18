import os
import numpy as np
import re
import matplotlib.pyplot as plt
import math

def get_current_folder():
    cur_dir1 = os.getcwd()  # works with Spyder
    cur_dir2 = cur_dir1.replace("\\", "/")
    return cur_dir1, cur_dir2

def main():
    epsilon = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
                        0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                        0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 
                        0.8, 0.85, 0.9])
    zlo = 0
    zhi = 200
    xlo = -75
    xhi = 75
    dr = 1

    NX = int((xhi-xlo)/dr)+2
    NZ = int((zhi-zlo)/dr)+2    
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
                        try:
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
                            
                            GZ[yg1][xg1] += w1x*w1y
                            GZ[yg2][xg1] += w1x*w2y
                            GZ[yg3][xg1] += w1x*w3y
                            
                            GZ[yg1][xg2] += w2x*w1y
                            GZ[yg2][xg2] += w2x*w2y
                            GZ[yg3][xg2] += w2x*w3y
                            
                            GZ[yg1][xg3] += w3x*w1y
                            GZ[yg2][xg3] += w3x*w2y
                            GZ[yg3][xg3] += w3x*w3y
                        except IndexError:
                            print("xi = %.2f, yi = %.2f, xg2 = %d, yg2 = %d" % (datax,datay,xg2,yg2) )
                        except:
                            print("Something else went wrong")
                            

        # ---plotting business
        plt.show()
        plt.contourf(GX, GY, GZ, 20, cmap='YlGnBu')
        plt.colorbar()
        plt.title("Epsilon = %.3f" % Temp)
        plt.axes().set_aspect('equal', 'datalim')
        plt.savefig(out_image, dpi=1200)

    print("\nProgram finished.")


if __name__ == "__main__":
    main()
