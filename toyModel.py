#! python

__author__ = 'shahar'
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from itertools import product
from scipy.optimize import fsolve, broyden1, broyden2, newton_krylov
from scipy.optimize.nonlin import NoConvergence
import sys
from optparse import OptionParser

def getM1(M, N):
    '''
    Creating the matrix represents
    Circhhuf's law. The external voltage equals the some of all voltages
    from source to drain
    :param M: Number of lines in the array
    :param N: Number of rows in the array
    :return: Matrix representation of the equation (so M1*V=Vext is the
    equation)
    '''
    M1 = np.zeros((M//2 + 1, M*N))
    for i in range(M//2 + 1):
        M1[i, 2*N*i:2*N*i + N] = 1
    return reduceArray(M1, M, N, axis=1)

def getM2(M,N):
    '''
    Creating the matrix represents
    Circhhuf's law. The sum of voltages in any
    closed loop is zero
    :param M: Number of lines in the array
    :param N: Number of rows in the array
    :return: Matrix representation of the equations (so M2*V=0 is the
    equation)
    '''
    M2 = np.zeros(((M//2)*(N-1), M*N))
    for line in range((M//2)*N):
        i = (line // (N-1)) * 2
        j = (line % (N-1))
        if i < M - 2:
            if 0 < j < N-1:
                M2[line, [i*N + j, (i+1)*N + j-1, (i+1)*N+j, (i+2)*N + j]] = \
                  [1, -1,  1, -1]
            elif j == 0:
                M2[line, [i*N + j, (i+1)*N+j, (i+2)*N + j]] = \
                  [1,  1, -1]
    return reduceArray(M2, M, N, axis=1)

def getM3NonLinear(M,N):
    '''
    Creating the matrix represents
    Circhhuf's law. The sum of entering and exiting current for each
    crossection is zero.
    This function is for the case where the resistors are non-Ohmic.
    :param M: Number of lines in the array
    :param N: Number of rows in the array
    :return: Matrix representation of the equations.
    The equations are M3 * I(V)=0
    '''
    M3 = np.zeros(((M//2 + 1)*(N - 1), M*N))
    for line in range((M//2 + 1)*(N - 1)):
        i = (line // (N-1)) * 2
        j = line % (N- 1)
        if 0 < i < M - 1:
            if j < N - 1:
                M3[line, [(i-1)*N + j, i*N + j, i*N+j+1, (i+1)*N + j]] = \
                    [1, 1,  -1, -1]
        elif i == M-1:
            if j < N - 1:
                M3[line, [(i-1)*N + j, i*N + j, i*N+j+1]] = \
                    [1, 1,  -1]
        elif i == 0:
            if j < N - 1:
                M3[line, [i*N + j, i*N+j+1, (i+1)*N + j]] = \
                    [1,  -1, -1]

    return reduceArray(M3, M, N , axis=1)

def getM3Linear(R,M,N):
    '''
    Creating the matrix represents
    Circhhuf's law. The sum of entering and exiting current for each
    crossection is zero.
    This function is for the case where the resistors Ohmic.
    :param R: Matrix of the resistance in the array
    :return: Matrix representation of the equations.
    The equations are M3 * V=0
    '''
    invR = expandArray(1/R, M,N).reshape((M,N))
    M3 = np.zeros(((M//2 + 1)*(N - 1), M*N))
    for line in range((M//2 + 1)*(N - 1)):
        i = (line // (N-1))*2
        j = line % (N-1)
        if 0 < i < M - 1:
            if j < N - 1:
                M3[line, [(i-1)*N + j, i*N + j, i*N+j+1, (i+1)*N + j]] = \
                    [invR[i-1, j], invR[i,j],  -invR[i,j+1], -invR[i+1,j]]
        elif i == M-1:
            if j < N - 1:
                M3[line, [(i-1)*N + j, i*N + j, i*N+j+1]] = \
                    [invR[i-1,j], invR[i,j],  -invR[i,j+1]]
        elif i == 0:
            if j < N - 1:
                M3[line, [i*N + j, i*N+j+1, (i+1)*N + j]] = \
                    [invR[i,j],  -invR[i,j+1], -invR[i+1,j]]
    return reduceArray(M3, M, N, axis=1)

def getDynamicPart(M,N,R,C,Vth,Vext,lastV,dt):
    '''
    Creating the matrix represents
    Circhhuf's law. The sum of entering and exiting current for each
    crossection is zero.
    This function is for the case where the resistors Ohmic.
    :param R: Matrix of the resistance in the array
    :return: Matrix representation of the equations.
    The equations are M3 * V=0
    '''
    Icoeff = np.zeros(lastV.shape)
    invR = 1/R
    Cdt = C/dt
    Icoeff[lastV >= Vth] = invR[lastV >= Vth]
    Icoeff[lastV < Vth] = Cdt[lastV < Vth]
    additional_b = np.zeros(lastV.shape)
    additional_b[lastV < Vth] = Icoeff[lastV < Vth]*lastV[lastV < Vth]
    Icoeff = expandArray(Icoeff, M, N).reshape((M, N))
    additional_b = expandArray(additional_b, M, N).reshape((M, N))
    M3 = np.zeros(((M//2 + 1)*(N - 1), M*N))
    b = getb(M,N,Vext)
    for line in range((M//2 + 1)*(N - 1)):
        b_line = -(M//2 + 1)*(N - 1) + line
        i = (line // (N-1)) * 2
        j = line % (N-1)
        if 0 < i < M - 1:
            if j < N - 1:
                M3[line, [(i-1)*N + j, i*N + j, i*N+j+1, (i+1)*N + j]] = \
                    [Icoeff[i-1, j], Icoeff[i,j],  -Icoeff[i,j+1], -Icoeff[i+1, j]]
                b[b_line] = additional_b[i-1, j] + additional_b[i,j] -\
                            additional_b[i,j+1] - additional_b[i+1, j]
        elif i == M-1:
            if j < N - 1:
                M3[line, [(i-1)*N + j, i*N + j, i*N+j+1]] = \
                    [Icoeff[i-1,j], Icoeff[i,j],  -Icoeff[i,j+1]]
                b[b_line] = additional_b[i-1, j] + additional_b[i,j] -\
                            additional_b[i, j+1]
        elif i == 0:
            if j < N - 1:
                M3[line, [i*N + j, i*N+j+1, (i+1)*N + j]] = \
                    [Icoeff[i,j],  -Icoeff[i,j+1], -Icoeff[i+1,j]]
                b[b_line] = additional_b[i, j] - additional_b[i,j+1] -\
                            additional_b[i+1,j]
            else:
                M3[line, [i*N + j, (i+1)*N + j]] = \
                    [Icoeff[i,j], -Icoeff[i+1,j]]
                b[b_line] = additional_b[i, j] - additional_b[i+1,j]
    return reduceArray(M3, M, N, axis=1), b

def getDynamicI(R, C, Vth, V, lastV, dt):
    I = np.zeros(V.shape)
    I[lastV >= Vth] = V[lastV >= Vth]/R[lastV >= Vth]
    I[lastV < Vth] = (V[lastV < Vth] - lastV[lastV < Vth]) * C[lastV < Vth] \
                     / dt
    return I
def getb(M,N, Vext):
    '''
    Returns the non homeogeneus part for the equations
    system for the ohmic case.
    :param R: Mtraix of the array resistances
    :param Vext: The external voltage applied
    :return: vector b so the equations are M*V = b
    '''
    b = np.zeros((M*N- M//2, 1))
    b[:M//2 + 1] = Vext
    return b

def getI(V, n, R0, Vth):
    '''
    Calculates I for the non Ohmic resistors.
    :param V: Matrix of the Voltages in the array
    :param n: parameter for the I(V) relation (when n-> infinity this case
    merges with the "linear" case
    :param R0: Resistance of resistor when V>Vth in the limit n-> infinity
    :param Vth: Matrix of threshold voltages
    :return:
    '''
    I = V**n / (R0*(Vth**(n-1) + V**(n-1)))
    I[V==0] = 0
    return I

def getIprime(V, n, R0, Vth):
    Vth = Vth.flatten()
    Iprime = (V**(n-1)/R0)*\
             ((n*(Vth)**(n-1) + V**n-1) / (Vth**(n-1) + V**(n-1))**2)
    Iprime[V==0] = 0
    return Iprime

def getEquations(M1, M2, M3, n, R0, Vth, Vext):
    '''
    Creating the equations system f for the non linear case.
    :param M: Number of rows in the array
    :param N: Number of coloumns in the array
    :param n: parameter for the non-Ohmic resistance
    :param R0: the resistance for resistor above threshold voltage
    :param Vth: threshold voltagees
    :param Vext: external voltage
    :return: f so the equations are f(V) = 0.
    '''

    def f(V):
        eq1 = M1.dot(V) - Vext
        eq2 = M2.dot(V)
        eq3 = M3.dot(getI(V, n, R0, Vth))
        return np.hstack((eq1, eq2, eq3))

    return f

def getLinearEquations(M,N, M1, M2, M3, Vext):
    '''
    Creating the equations system for the linear case.
    :param R: Matrix of array resistances
    :param Vext: external voltage
    :return: M, b so the equations are M*V = b
    '''
    Mat = np.vstack((M1, M2, M3))
    b = getb(M,N, Vext)
    return Mat, b

def calcDynamicVoltages(M, N, M1, M2, R, C, Vth, Vext, lastV, dt, repeat=1,
                        voltageTime = 1):
    lastV = lastV.flatten()
    Ilist = []
    for time in range(repeat*voltageTime):
        M3,b = getDynamicPart(M, N, R, C, Vth, Vext, lastV, dt)
        Mat = np.vstack((M1, M2, M3))
        V = np.linalg.solve(Mat, b)
        V = V.flatten()
        if time >= (voltageTime-1)*repeat:
            Ilist.append(getDynamicI(R, C, Vth, V, lastV, dt))
        lastV = V
    return Ilist, V

def calcLinearCurrents(M,N, M1, M2, R, R0, Rinf, Vth, Vext, repeat):
    '''
    Calculating voltages for the Ohmic case.
    :param R: Resistances matrix
    :param R0: resistance of resistor above threshold
    :param Rinf: resistance of resistor below threshold
    :param Vth: threshold voltages
    :param Vext: external voltage
    :return: V,R. Voltages and Resistances.
    '''
    print('Start Calculating')
    M3 = getM3Linear(R,M,N)
    Mat, b = getLinearEquations(M,N, M1, M2, M3, Vext)
    alreadyChanged = np.zeros(R.shape)
    V = np.linalg.solve(Mat, b).flatten()
    Rnew = updateR(R, V, Vth, R0, Rinf)
    # if Rnew is different from are the solution we found is not consistent
    #  with the threshold rules
    res = []
    while (Rnew != R).any():
        alreadyChanged = alreadyChanged + (Rnew != R).astype(np.int)
        R = Rnew
        M3 = getM3Linear(R,M,N)
        Mat, b = getLinearEquations(M,N,M1, M2,M3 ,Vext)
        V = np.linalg.solve(Mat, b).flatten()
        Rnew = updateR(R, V, Vth, R0, Rinf)
        # if np.max(alreadyChanged) > 10:
        #     res.append(V/R)
        #     if len(res) > repeat:
        #         break
    if len(res) == 0:
        res.append(V/R)
    return res, Rnew

def getJacobian(M1, M2, M3,n,R0,Vth):
    def fprime(V):
        nonLinearPart = M3.dot(np.diag(getIprime(V, n, R0, Vth)))
        if M2.size:
            return np.vstack((M1, M2, nonLinearPart))
        else: # in the case of one dimensional array
            return np.vstack((M1, nonLinearPart))
    return fprime

def calcNonLinearVoltages(M, N,M1, M2,M3, n, R0, Vth, Vext, Vini):
    '''
    Calculating Voltages for the non-Ohmic case.
    :param M: Number of rows in the array
    :param N: Number of coloumns in the array
    :param n: parameter for the non-Ohmic resistance
    :param R0: the resistance for resistor above threshold voltage
    :param Vth: threshold voltagees
    :param Vext: external voltage
    :param Vini: Initial guess for voltages
    :return: Voltages of the array
    '''
    f = getEquations(M1, M2, M3, n, R0, Vth.flatten(), Vext)
    # fprime = getJacobian(M1, M2, M3, n, R0, Vth.flatten())
    # V, info, flag, msg = fsolve(f, Vini, fprime=fprime, full_output=True)
    # if flag !=1:
    #     print(msg)
    #     raise NoConvergence
    V = broyden1(f, Vini)
    return V


def checkV(R, V, Vth, R0, Rinf):
    '''
    Checks if V found is consistent with the thresholds.
    Used for the Ohmic case
    :param V: Voltages
    :param Vth: Threshold voltages
    :param R0: Resistance of resistor above threshold
    :param Rinf: Resistance of resistor below threshold
    :return: True if the voltages are consistent with threshold
    '''
    above = np.abs(V) >= Vth
    below = np.abs(V) < Vth
    falseConducting = np.logical_and(below, R == R0)
    falseInsulating = np.logical_and(above, R == Rinf)
    return not (falseConducting.any() or falseInsulating.any())

def generateAllR(R0, Rinf, M, N):
    '''
    Generates all possible R combinations
    :param R0: Resistance of resistor above threshold
    :param Rinf: Resistance of resistor below threshold
    :param M: Number of rows in the array
    :param N: Number of coloumns in the array
    :return:
    '''
    for R in product([R0, Rinf], repeat=N*M):
        yield np.reshape(np.array(R), (M,N))

def updateR(R, V, Vth, R0, Rinf):
    '''
    Update R matrix so it will match threshold conditions
    Workes by randomly choosing resistor wich are not consistent and
    change them.
    :param R: Current resistances
    :param V: Current Voltages
    :param Vth: Threshold voltages
    :param R0: Resistance of resistor above threshold
    :param Rinf: Resistance of resistor below threshold
    :return: The new matrix. If all resistances are consistent R will be
    returned.
    '''
    V = V.flatten()
    Rnew = np.copy(R)
    above = np.abs(V) >= Vth
    below = np.abs(V) < Vth
    falseConducting = np.array(np.where(np.logical_and(below, R == R0)))
    falseInsulating = np.array(np.where(np.logical_and(above, R == Rinf)))
    falseConducting = falseConducting[0]
    falseInsulating = falseInsulating[0]
    if falseConducting.size:
        turnToInsulator = np.random.randint(0, falseConducting.size)
        Rnew[falseConducting] = Rinf
    if falseInsulating.size:
        turnToConductor = np.random.randint(0, falseInsulating.size)
        Rnew[falseInsulating] = R0[falseInsulating]
    print("false conducting: " + str(falseConducting.size))
    print("false insulating: " + str(falseInsulating.size))
    return Rnew

def reduceArray(a, M, N, axis=0):
    return np.delete(a,np.arange(1,M,2)*N + N-1, axis=axis)

def expandArray(a, M , N, axis=0):
    return np.insert(a,np.arange(2*N-1, (2*N-1)*M//2, 2*N-1),[0], axis=axis)

def calcIVcurve(M, N, n, Rinf, Vth, Vmin, Vmax, stepNum, R,  C, dt, repeat,
                method="nonlinear", voltageTime=1):
    '''
    Calculates the relation between external current and voltage
    :param M: Number of rows in the array
    :param N: Number of coloumns in the array
    :param n: parameter for the non-Ohmic resistance
    :param R0: the resistance for resistor above threshold voltage
    :param Rinf: tM, N, R0, Rinf, Vth, Vmaxhe resistance for resistor below threshold voltage
    :param Vth: threshold voltages
    :param Vmax: Maximum external voltage
    :param linear: If Ture uses the Ohmic version. The non-linear version
    if false
    :param n: Parameter for the non-Ohmic version
    :return: Iext, Vext, ims list of external voltage and currents and
    images for animation.
    '''
    M1 = getM1(M,N)
    M2 = getM2(M,N)
    M3nonLinear = getM3NonLinear(M,N)
    results = []
    Iext = []
    Vext = np.linspace(Vmin,Vmax, num=stepNum)
    Vext = np.hstack((Vext, np.flip(Vext, axis=0)))
    R = reduceArray(R.flatten(), M, N)
    C = reduceArray(C.flatten(), M, N)
    currR = np.ones((M*N-M//2,))*Rinf
    V = np.zeros((M*N,))
    V = np.delete(V,np.arange(1,M,2)*N + N-1)
    for index, currV in enumerate(Vext):
        try:
            if method=="linear":
                Ilist, newR= calcLinearCurrents(M,N, M1, M2, currR,R,
                                               Rinf, Vth, currV, repeat)
                currR = newR
            elif method == "dynamic":
                Ilist, V = calcDynamicVoltages(M, N, M1, M2, R, C, Vth,
                                               currV, V, dt, repeat,
                                               voltageTime=voltageTime)
            elif method == "nonlinear":
                V = calcNonLinearVoltages(M, N,M1, M2,M3nonLinear, n, R,
                                          Vth, currV, V)
                Ilist = [getI(V,n,R, Vth)]
            else:
                print("no such method: " + method)
                exit(1)
        except NoConvergence:
            results.append((None,None))
            Iext.append(-999)
            print("failed for V=" + str(currV))
            continue
        results.extend([(I, currV) for I in Ilist])
        Iext.append(np.mean([np.sum(I[(N-1)::(2*N-1)]) for I in Ilist]))
    Iext = np.array(Iext)
    Vext = np.array(Vext)
    return Iext[Iext != -999], Vext[Iext != -999], results

def plotCurrent(results, Vmax, M,N):
    '''
    updating the plot to current currents map
    :return: image for animation
    '''
    ims = []
    num = 0
    for result in results:
        I, Vext = result
        if I is None:
            continue
        I = expandArray(I,M,N)
        I = I.reshape((M,N))
        J = np.zeros(((M//2)*3+1,N*3))
        for i in range(2):
            J[0:(M//2)*3+1:3, i:3*N:3] = I[0:M:2,:]
        for i in range(2):
            J[i+1:(M//2)*3+1:3, 2:3*N:3] = I[1:M:2,:]
        im = plt.imshow(J, vmin=-Vmax/(M*M*N), vmax=Vmax/(M*M*N),
        cmap='plasma', aspect='equal',animated=True)
        tx = plt.annotate('Vext = ' + str(Vext), (1,1))
        ims.append([im, tx])
        num += 1
    print(str(num))
    return ims

def getDisorderedValues(M, N, mean, std, nozero=False, distribution="normal"):
    minval = 0.01 if nozero else 0
    if distribution == "normal":
        arr = np.clip(np.random.normal(loc=mean, scale=std, size=(M,N)), minval,
                                       mean+3*std)
    elif distribution == "exponential":
        arr = np.random.exponential(scale=std, size = (M,N)) + (mean - std)
    else:
        print("The distribution " + distribution + " is not implemented yet")
        raise NotImplementedError
    return arr

def plotVth(Vth):
    M,N = Vth.shape
    copyVth = np.copy(Vth)
    copyVth[1:M:2, N-1] = 0
    J = np.zeros(((M//2)*3+1,N*3))
    for i in range(2):
        J[0:(M//2)*3+1:3, i:3*N:3] = copyVth[0:M:2,:]
    for i in range(2):
        J[i+1:(M//2)*3+1:3, 2:3*N:3] = copyVth[1:M:2,:]
    J[J < 0.001] = 0
    plt.imshow(J, vmin=0, vmax=np.max(copyVth), cmap='jet', aspect='auto')
    plt.colorbar()

def runModel(N,M,n,Vmin, Vmax, numSteps, VthMean,
             VthStd, RMean, RStd, CMean, CStd, dt,  method,
             fileNamesPostfix, repeat, distribution="normal",plot = False,
             voltageTime=1):
    if RStd == 0:
        R = np.ones((M, N))*RMean
    else:
        R = getDisorderedValues(M, N, RMean, RStd, nozero=True)
    if CStd == 0:
        C = np.ones((M, N))*CMean
    else:
        C = getDisorderedValues(M, N, CMean, CStd, nozero=True)
    Rinf = np.max(R) * 1000
    Vth = getDisorderedValues(M, N, VthMean, VthStd,
                              distribution=distribution)
    fig1 = plt.figure()
    plotVth(Vth)
    plt.title('Trheshold voltages for orderParam = ' + str(VthStd))
    plt.savefig('ThresholdVoltages_orderParam_' + str(VthStd) + "_vtime" +
                str(voltageTime) + "_size=" + str(N) + "X" + str(M) +
                fileNamesPostfix +".png")
    plt.close(fig1)
    Vth = reduceArray(Vth.flatten(), M, N)
    Iext, Vext, results = calcIVcurve(M, N, n, Rinf, Vth, Vmin, Vmax,
                                      numSteps, R,  C, dt,
                                      repeat, method=method,
                                      voltageTime=voltageTime)

    if plot:
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=24, bitrate=1800)
        fig2 = plt.figure()
        plt.imshow(np.zeros((M,N)), vmin=-Vmax/np.min(R), vmax=Vmax/np.min(R),
                   animated=True, cmap='plasma', aspect='equal')
        plt.colorbar()
        ims = plotCurrent(results, Vmax, M,N)
        im_ani = animation.ArtistAnimation(fig2, ims, interval=100,
                                           repeat_delay=1000, blit=True)
        im_ani.save('CurrentsAnimation_orderParam_' + str(VthStd) + "_vtime" +
                str(voltageTime) + "_size=" + str(N) + "X" + str(M) +
                fileNamesPostfix + '.mp4',
                writer=writer)
        plt.close(fig2)
    fig3 = plt.figure()
    x = len(Vext)
    plt.plot(Vext[:x//2], Iext[:x//2],
             'b.',Vext[x//2:], Iext[x//2:], 'r.')
    plt.xlabel('Voltage')
    plt.ylabel('Current')
    plt.title('IV for orderParam = ' + str(VthStd))
    plt.savefig('IV_orderParam_' +  str(VthStd) +"_vtime" +
                str(voltageTime) + "_size=" + str(N) + "X" + str(M) +
                fileNamesPostfix +".png")
    plt.close(fig3)
    with open('runningParameters_orderParam_' +  str(VthStd) + "_vtime" +
                str(voltageTime) + "_size=" + str(N) + "X" + str(M) +
                fileNamesPostfix +".txt",mode='w') as f:
        f.write("N: " + str(N) + "\nM: " + str((M+1)//2)+ "\nn: " +
                str(n) + "\nVmin: " + str(Vmin) +  "\nVmax: " +
                str(Vmax) + "\n step num: " + str(numSteps) + "\nVth mean: "+
                str(VthMean) + "\n Vth std: " + str(VthStd) +
                "\n R mean: " + str(RMean)+"\n R std: " + str(RStd)+
                "\nC mean: " + str(CMean) + "\n C std: " + str(CStd) +
                "\ndt: " + str(dt) + "\nmethod: " + method + "\nrepesats: " +
                str(repeat))

def getOptions():
    parser = OptionParser(usage= "usage: %prog [options]")

    parser.add_option("-M", "--height", dest="M", help="number of lines in "
                      "the array [default: %default]", default=1, type=int)
    parser.add_option("-N", "--width", dest="N", help="number of columns in "
                      "the array [default: %default]", default=1, type=int)
    parser.add_option("--vmin", dest="Vmin", help="minimum external voltage"
                      " [default: %default]", default=0, type=float)
    parser.add_option("--vmax", dest="Vmax", help="maximum external voltage"
                      " [default: %default]", default=100, type=float)
    parser.add_option("--stepnum", dest="stepNum", help="number of voltage "
                      "steps[default: %default]", default=1000, type=int)

    parser.add_option("--vthmean", dest="VthMean", help="average threshold "
                      "voltage value[default: %default]", default=100,
                      type=float)
    parser.add_option("--vthstd", dest="VthStd", help="standart deviation "
                      "threshold voltages values[default: %default]",
                      default=[], type=float, action="append")
    parser.add_option("--vtimes", dest="voltageTimes",
                      help="how many times to run dt befotr changing voltage"
                           " [default: %default]",
                      default=[], type=int, action="append")
    parser.add_option("--rmean", dest="RMean", help="mean "
                      "resistance values[default: %default]",
                      default=1, type=float)
    parser.add_option("--rstd", dest="RStd", help="standart deviation "
                      "of resistances values[default: %default]",
                      default=0, type=float)
    parser.add_option("--cmean", dest="CMean", help="mean "
                      "capacitance values[default: %default]",
                      default=1, type=float)
    parser.add_option("--cstd", dest="CStd", help="standart deviation "
                      "of capacitance values[default: %default]",
                      default=0, type=float)
    parser.add_option("-n", dest="n", help="n value used in nonlinear "
                      "resistance [default:%default]",
                      default=5, type=float)
    parser.add_option("--dt", dest="dt", help="time step size"
                      " [default:%default]", default=0, type=float)
    parser.add_option("--method", dest="method", help="method for "
                      "calculation, from (linear, nonlinear, dynamic) "
                      "[default: %default]", default="linear")
    parser.add_option("--distribution", dest="distribution",
                      help="probability distribution from wich the Vth "
                           "values are drawn [default: %default]",
                      default="normal")
    parser.add_option("--postfix", dest="fileNamesPostfix", help="optional "
                      "addition to output filenames", default='')
    parser.add_option("--repeatnum", dest="repeat", help="number of times "
                      "to repeat each voltage for averaging (for linear "
                      "meathod) [default: %default]", default=1, type=int)
    parser.add_option("--plotcurrent", dest="plot", help="if true a clip "
                      "with currnets will be plotted [Default: %default]",
                      default=False, action='store_true')
    return parser.parse_args()

if __name__ == "__main__":
    options, args = getOptions()
    for vthStd in options.VthStd:
        for voltageTime in options.voltageTimes:
            runModel(options.N, options.M*2-1, options.n,
                options.Vmin, options.Vmax, options.stepNum, options.VthMean,
                vthStd, options.RMean, options.RStd, options.CMean,
                options.CStd, options.dt,  options.method,
                options.fileNamesPostfix, options.repeat,
                 distribution=options.distribution, plot=options.plot,
                 voltageTime=voltageTime)
