print('Program started')
import os
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


with open(os.path.join(__location__, 'parameters.txt')) as fileID:
    for line in fileID:
        temp = line.split();
        if len(temp) > 0:
            if temp[0] == "Current":
                print("%s %s" % (temp[3],temp[6]))
        
print('Program finished')
