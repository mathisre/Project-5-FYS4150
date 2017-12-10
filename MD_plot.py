from sys import argv
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["savefig.directory"] = "/home/mathisre/Dropbox/Uni/Images/Proj5"

if len(argv) != 5:
    print "This program uses command line arguments to read files titled temp_diff_x.txt."
    print "Write location of files, start x, final x and increments."
    exit()


n = 100000 #just a large enough number
tempInit = int(argv[2])
tempFinal= int(argv[3])
tempStep = int(argv[4])


DiffC = []
DiffT_rel = []
DiffT = []


for i in range(tempInit, tempFinal+tempStep, tempStep):
    file = open(argv[1] + 'temp_diff_' + str(i) + '.txt', 'r')
    with file as filename:
        lines = [line.split() for line in filename]

    temperature = np.zeros(n)
    time = np.zeros(n)
    diffConstant = np.zeros(n)

    for k in range(len(lines)): #read data into array
        time[k]         = float(lines[k][0])
        temperature[k]  = float(lines[k][1])
        diffConstant[k] = float(lines[k][2])
    file.close()

    #Remove empty entries from arrays
    temperature = np.trim_zeros(temperature , 'b')
    rel_temperature = temperature / temperature[0]
    time= np.trim_zeros(time,'b')
    diffConstant = np.trim_zeros(diffConstant,'b')

    marker = "Initial temperature = %02d K" % (np.floor(temperature[0]))
    plt.figure(0)
    plt.plot(time, rel_temperature, label=marker)

    plt.figure(1)
    plt.plot(time[1:-1],diffConstant[1:-1], label = marker)


    #Create mean values of equilibrium temperature and diffusion constant
    N = 1000
    meanDiffC = 0
    meanDiffT = 0
    for j in range(1, N + 1):
        meanDiffC += diffConstant[-j]
        meanDiffT += temperature[-j]
    meanDiffC /= N
    meanDiffT /= N

    DiffC.append(meanDiffC)
    DiffT.append(meanDiffT)
    DiffT_rel.append(meanDiffT/temperature[0])

plt.figure(0)
plt.xlabel('Time in units of 10^-13 seconds')
plt.ylabel('$T/T_i$')
plt.title('Temperature time evolution')
#plt.legend(loc = 'best')
plt.grid(True)
#plt.axis([0, 30, 0,1])

plt.figure(1)
plt.xlabel('Time in units of 10^-13 seconds')
plt.ylabel('Diffusion constant [$m^2/s$]')
plt.title('Diffusion constant time evolution')
#plt.legend(loc = 'best')
plt.grid(True)


plt.figure(2)
plt.scatter(DiffT, DiffC,marker = 'x')
plt.xlabel('Temperature [K]')
plt.ylabel('Diffusion constant [$m^2/s$]')
plt.title('Equilibrium diffusion constant')
plt.axis([50, 450, 0, 6*10**-11])
plt.grid(True)


plt.figure(3)
plt.scatter(DiffT , DiffT_rel,marker ='x')
plt.xlabel('Temperature [K]')
plt.ylabel('$T/T_i$')
plt.title('Variance of equilibrium temperature ratio')
plt.grid(True)
plt.axis([50, 450, 0.44, 0.55])

plt.show()