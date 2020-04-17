import os
import numpy as np
import re
import matplotlib.pyplot as plt


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


def main():
    epsilon = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
                        0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    # zlo = 0;    zhi = 200;    xlo = -75;    xhi = 75;

    dir1, dir2 = get_current_folder()

    for Temp in epsilon:
        in_file = dir2 + "/data/dump_" + str(Temp) + ".xyz"
        out_file = dir2 + "/results/result_" + str(Temp) + ".txt"
        out_image = dir2 + "/results/image_" + str(Temp) + ".png"
        print("File to read is \n", in_file)
        print("File to write is \n", out_file)

        zlo = 0
        zhi = 200
        xlo = -75
        xhi = 75
        dr = 1
        NX = int((xhi-xlo)/dr)
        NZ = int((zhi-zlo)/dr)
        res = np.zeros([NZ, NX])  # note that columns are x and rows are z

        print("\nPlease wait while I am working...")

        with open(in_file) as file_in:
            for line in file_in:
                temp1 = re.findall(
                    r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)
                # An expected line should be 2 -47.25 -74.25 0 [type x y z]
                if len(temp1) == 4:  # for a valid entry
                    if float(temp1[0]) == 1:  # for type 1, water beads
                        locx = int((float(temp1[1])-xlo)/dr)
                        locz = int((float(temp1[3])-zlo)/dr)
                        res[locz][locx] += 1
        # ---write the optional data to file
        with open(out_file, "w") as file_out:
            for i in range(NZ):
                for j in range(NX):
                    file_out.write("%.1f " % res[i][j])
                file_out.write("\n")

        maxval = np.amax(res)
        print("Max value: %f" % maxval)

        # ---plotting business
        x, y = np.meshgrid(np.linspace(xlo, xhi, NX),
                           np.linspace(zlo, zhi, NZ))
        plt.show()
        #plt.contourf(x, y, res, 10, cmap='binary')
        plt.imshow(res, extent=[xlo, xhi, zlo, zhi],
                   origin='lower', cmap='jet')
        plt.colorbar()
        plt.axis(aspect='image')
        plt.savefig(out_image, dpi=1200)

    print("\nProgram finished.")


if __name__ == "__main__":
    main()
