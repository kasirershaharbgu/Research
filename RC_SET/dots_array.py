__author__ = 'shahar'

import numpy as np
from matplotlib import pyplot as plt
from optparse import OptionParser
import os
from ast import literal_eval

def flattenToColumn(a):
    return a.reshape((a.size,1))

class NoElectronsOnDot(RuntimeError):
    pass

class DotArray:
    """
    An array of quantum dots connected to external voltage from left to
    right and to gate voltage
    """

    def __init__(self, rows, columns, Vext, VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv):
        """
        Creates new array of quantum dots
        :param rows: number of rows (int)
        :param columns: number of columns (int)
        :param Vext: external voltage (double)
        :param Q0: np array of initial charges (rowsXcolumns double array)
        :param n0: np array of initial electrons (rowsXcolumns double array)
        :param VG: voltage of the gates (rowsXcolumns double array)
        :param Ch: horizotal capacitances ((rowsXcolumns+1 double array)
        :param Cv: vertical capacitances (rows-1Xcolumns double array)
        :param Rh: vertical tunnelling ressistances (rowsXcolumns+1 double array)
        :param Rv: horizontal tunnelling ressistances (rowsXcolumns+1 double array)
        :param RG: ground resistances (rowsXcolumns double array)
        """
        self.rows = rows
        self.columns = columns
        self.VL = Vext/2
        self.VR = -Vext/2
        self.VG = VG
        self.Q = Q0
        self.n = n0
        self.CG = CG
        self.RG = RG
        self.Ch = Ch
        self.Cv = np.vstack((np.zeros((1,columns)),Cv))
        self.Rh = Rh
        self.Rv = Rv
        self.totalChargePassed = 0
        self.createCapacitanceMatrix()
        self.setDiagonalizedJ()
        self.setConstWork()

    def getRows(self):
        return self.rows

    def getColumns(self):
        return self.columns

    def resetCharge(self):
        """
        Reset the counter for total charge passed
        """
        self.totalChargePassed = 0
        return True

    def changeVext(self, newVext):
        self.VL = newVext/2
        self.VR = -newVext/2
        return True

    def getCharge(self):
        return self.totalChargePassed

    def createCapacitanceMatrix(self):
        """
        Creates the inverse capacitance matrix
        """
        diagonal = self.Ch[:,:-1] + self.Ch[:,1:] + self.Cv + np.roll(self.Cv, -1, axis=0)
        second_diagonal = np.copy(self.Ch[:,1:])
        second_diagonal[:,-1] = 0
        second_diagonal = second_diagonal.flatten()
        second_diagonal = second_diagonal[:-1]
        n_diagonal = np.copy(self.Cv[1:,:])
        C_mat = np.diagflat(diagonal) - np.diagflat(second_diagonal,k=1) - np.diagflat(second_diagonal,k=-1)\
                -np.diagflat(n_diagonal, k=self.columns) -  np.diagflat(n_diagonal, k=-self.columns)
        self.invC = np.linalg.inv(C_mat)
        return True

    def getNprime(self):
        left_part = np.copy(self.Ch)
        left_part[:,1:] = 0
        right_part = np.copy(self.Ch)
        right_part[:,:-1] = 0
        return flattenToColumn(self.n + left_part*self.VL + right_part*self.VR)

    def getbVector(self):
        res = -self.invC.dot(self.getNprime()) + flattenToColumn(self.VG)
        return res / np.repeat(flattenToColumn(self.RG),res.shape[1],axis=1)

    def getJmatrix(self):
        res = self.invC + np.diagflat(1/self.CG)
        return res / np.repeat(flattenToColumn(self.RG),res.shape[1],axis=1)

    def setDiagonalizedJ(self):
        self._JeigenValues, self._JeigenVectors = np.linalg.eig(self.getJmatrix())
        self._JeigenVectorsInv = np.linalg.inv(self.JeigenVectors)
        return True

    def developeQ(self, dt):
        Q0 = self._JeigenVectorsInv.dot(self.Q.T)
        b = self.getbVector()
        exponent = np.exp(self._JeigenValues*dt)
        self.Q = self._JeigenVectors.dot(Q0*exponent +(
             b/self._JeigenValues)*(exponent - 1)).T
        return True

    def setConstWork(self):
        invCDiagMat = np.diag(self.invC).reshape((self.rows,self.columns))
        lowerCDiag = np.pad(np.diag(self.invC,k=-1),((0,1),)).reshape((self.rows,self.columns))
        lowerCDiag[:,:-1] = 0
        lowerCDiag  = np.pad(lowerCDiag,((1,0),(0,0)))
        upperCDiag = np.pad(np.diag(self.invC,k=1),((0,1),)).reshape((self.rows,self.columns))
        upperCDiag[:,:-1] = 0
        upperCDiag  = np.pad(upperCDiag,((1,0),(0,0)))
        commonHorz = np.pad(invCDiagMat,((0,0),(1,0))) + np.pad(invCDiagMat,((0,0),(0,1)))\
                        - lowerCDiag - upperCDiag

        lowerNCDiag = np.diag(self.invC,k=-self.columns).reshape((self.rows-1,self.columns))
        upperNCDiag = np.diag(self.invC,k=self.columns).reshape((self.rows-1,self.columns))
        commonVert = invCDiagMat[1:,:] + invCDiagMat[:-1,:] - lowerNCDiag - upperNCDiag

        additionalLeft = np.zeros((self.rows,self.columns+1))
        additionalLeft[:,0] = self.VL
        additionalLeft[:,-1] = -self.VR
        self.leftConstWork = (0.5*commonHorz + additionalLeft).flatten()

        additionalRight = -additionalLeft
        self.rightConstWork = (0.5*commonHorz + additionalRight).flatten()

        self.upConstWork = (0.5*commonVert).flatten()
        self.downConstWork = np.copy(self.upConstWork)
        return True

    def getWork(self):
        q = self.getNprime() + flattenToColumn(self.QG)
        firstHorzMat = np.zeros(((self.columns+1)*self.rows,self.columns*self.rows))
        secondHorzMat = np.zeros(firstHorzMat.shape)
        firstLocations = np.ones(((self.columns+1)*self.rows,))
        secondLocations = np.ones(firstLocations.shape)
        firstLocations[self.columns:self.columns+1:] = 0
        secondLocations[0:self.columns+1:] = 0
        firstHorzMat[np.astype(firstLocations,np.bool),:] = self.invC
        secondHorzMat[np.astype(secondLocations,np.bool),:] = self.invC
        variableRightWork = ((firstHorzMat - secondHorzMat).dot(q)).flatten()
        variableLeftWork = -variableRightWork

        firstVertMat = self.invC[:-self.columns,:]
        secondVertMat = self.invC[self.columns:,:]
        variableUpWork = ((firstVertMat - secondVertMat).dot(q)).flatten()
        variableDownWork = -variableUpWork

        return variableRightWork + self.rightConstWork, variableLeftWork + self.leftConstWork,\
               variableDownWork + self.downConstWork, variableUpWork + self.upConstWork

    def getRates(self):
        rightWork, leftWork, downWork, upWork = self.getWork()
        horzResistance = self.Rh.flatten()
        vertResistance = self.Rv.flatten()
        rightWork[rightWork > 0] = 0
        leftWork[leftWork > 0] = 0
        downWork[downWork > 0] = 0
        upWork[upWork > 0] = 0
        return -rightWork/horzResistance, -leftWork/horzResistance,\
               -downWork/vertResistance,-upWork/vertResistance


    def getprobabilities(self, dt):
        rates = self.getRates()
        ratesCombined = np.hstack(rates)
        prob = 1 - np.exp(-ratesCombined*dt)
        return prob

    def tunnel(self, fromDot, toDot):
        if self.columns > fromDot[1] >= 0:
            self.n[fromDot] -=1
        if self.columns > toDot[1] >= 0:
            self.n[toDot] += 1
        if toDot[1] == self.columns:
            self.totalChargePassed += 1
        elif fromDot[1] == self.columns:
            self.totalChargePassed -= 1
        return True

    def printState(self):
        for i in range(self.rows):
            for j in range(self.columns):
                print("At dot (" + str(i) + ", " + str(j) + ") the state "
                      "is: n= " + str(self.n[i,j]) + " and QG = " + str(
                    self.Q[i,j]))

    def executeStep(self, ind):
            horzSize = (self.rows*(self.columns+1))**2
            vertSize = ((self.rows - 1)*self.columns)**2
            if ind <  horzSize: # tunnel right:
                fromDot = (ind//(self.columns+1), ind%(self.columns+1)-1)
                toDot = fromDot + (0,1)
            elif ind < horzSize*2:
                ind -= horzSize
                fromDot = (ind//(self.columns+1), ind%(self.columns+1))
                toDot = fromDot + (0,-1)
            elif ind < horzSize*2 + vertSize:
                ind -= horzSize*2
                fromDot = (ind//self.columns, ind%self.columns)
                toDot = fromDot + (1,0)
            else:
                ind -= (horzSize*2 + vertSize)
                fromDot = (ind//self.columns+1, ind%self.columns)
                toDot = fromDot + (-1,0)
            self.tunnel(fromDot, toDot)

MINPROBLOG = np.log(1/0.7)
class Simulator:
    """
    Preforms Markov-Chain Monte-Carlo simulation on an array of quantum dots
    :param rows: number of rows in the array
    :param columns: number of columns in the array
    """
    def __init__(self, rows, columns, Vext0, VG, dt,
                 Q0, n0, CG, RG, Ch, Cv, Rh, Rv):
        self.dotArray = DotArray(rows, columns, Vext0, VG, Q0, n0, CG, RG,
                                 Ch, Cv, Rh, Rv)
        self.t = 0
        self.dt = dt
        self.Vext = Vext0
        # for debug
    #     self.randomGen = self.getRandom()
    #
    # def getRandom(self):
    #     numbers = [0.01,0.001,0.001,0.76,0.002]
    #     idx = -1
    #     while True:
    #         idx += 1
    #         if idx == len(numbers):
    #             idx = 0
    #         yield numbers[idx]

    def executeStep(self, printState=False):
        r = np.random.ranf()
        # for debug
        # r = next(self.randomGen)
        probNotOk = True
        while (probNotOk):
            self.dotArray.developeQ(self.dt)
            prob = self.dotArray.getprobabilities(self.dt)
            cumProb = np.cumsum(prob)
            if cumProb[-1] > 1:
                print("Warning: total prob > 1 needed smaller dt, reducing dt")
                self.dt /= 10
            else:
                probNotOk = False
        actionInd = np.searchsorted(cumProb, r)
        if actionInd < prob.size:
            self.dotArray.executeStep(actionInd)
        if printState:
            self.printState()

    def calcCurrent(self, t,print=False):
        # if fullOutput:
        #     n = []
        #     Q = []
        #TODO: find smart way to know we have reach steady state
        steadySteps = 1000

        # Reaching equilibrium before calculating current
        for steps in steadySteps:
            self.executeStep(printState=print)
            self.dotArray.resetCharge()
        # Now when we are in steady state actually measuring current
        for steps in range(int(t // self.dt)):
            self.executeStep(printState=print)
        current = self.dotArray.getCharge() / t
        return current

    def calcIV(self, Vmax, Vstep, tStep, print=False, fullOutput=False):
        V = []
        I = []

        V0 = self.Vext
        while(self.Vext < Vmax):
            if fullOutput:
                current, stepn, stepQ = self.calcCurrent(tStep,print=print)
                I.append(current)
            else:
                I.append(self.calcCurrent(tStep, print=print))
            V.append(self.Vext)
            self.Vext += Vstep
            self.t += tStep
            self.dotArray.changeVext(self.Vext)
            self.dotArray.resetCharge()
        while(self.Vext > V0):
            if fullOutput:
                current, stepn, stepQ = self.calcCurrent(tStep,print=print)
                I.append(current)
            else:
                I.append(self.calcCurrent(tStep, print=print))
            V.append(self.Vext)
            self.Vext -= Vstep
            self.t += tStep
            self.dotArray.changeVext(self.Vext)
            self.dotArray.resetCharge()
        return np.array(I), np.array(V)

    def printState(self):
        print("At t= " + str(self.t) + " Vext = " + str(self.Vext) + ":")
        self.dotArray.printState()


def run1DSimulation(Vext0, VG0, dt, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, rows, columns,
                    Vmax, Vstep, tStep,repeats=1, savePath=".",
                    fileName="", fullOutput=False, printState=False):

    basePath = os.path.join(savePath, fileName)
    Is = []

    for repeat in range(repeats):
        simulator = Simulator(rows, columns, Vext0, np.array(VG0), dt,
                              np.array(Q0), np.array(n0), np.array(CG),
                              np.array(RG), np.array(Ch), np.array(Cv),
                              np.array(Rh), np.array(Rv))
        out = simulator.calcIV(Vmax, Vstep, tStep,
                                fullOutput=fullOutput, print=printState)
        I = out[0]
        V = out[1]
        # if fullOutput:
        #     n = out[2]
        #     Q = out[3]
        Is.append(I)
    avgI = np.mean(np.array(Is), axis=0)
    # if fullOutput:
    #     fig = plt.figure()
    #     plt.plot(n, '.')
    #     plt.savefig(basePath + "_n.png")
    #     np.save(basePath + "_n.bin", n)
    #     plt.close(fig)
    #
    #     fig = plt.figure()
    #     plt.plot(Q, '.')
    #     plt.savefig(basePath + "_Q.png")
    #     np.save(basePath + "_Q.bin", Q)
    #     plt.close(fig)

    fig = plt.figure()
    plt.plot(V[:V.size//2], avgI[:I.size//2], '.b',V[V.size//2:],
             avgI[I.size//2:],
             '.r')

    plt.savefig(basePath + "_IV.png")
    np.save(basePath + "_I.bin", avgI)
    np.save(basePath + "_V.bin", V)
    plt.close(fig)
    # if fullOutput:
    #     return I, V, n, Q
    return I, V

def getOptions():
    parser = OptionParser(usage= "usage: %prog [options]")

    # Normal parameters
    parser.add_option("-M", "--height", dest="M", help="number of lines in "
                      "the array [default: %default]", default=1, type=int)
    parser.add_option("-N", "--width", dest="N", help="number of columns in "
                      "the array [default: %default]", default=1, type=int)
    parser.add_option("--vmin", dest="Vmin", help="minimum external voltage"
                      " [default: %default]", default=0, type=float)
    parser.add_option("--vmax", dest="Vmax", help="maximum external voltage"
                      " [default: %default]", default=10, type=float)
    parser.add_option("--vstep", dest="vStep", help="size of voltage step["
                    "default: %default]", default=1, type=float)
    parser.add_option("--tstep", dest="tStep", help="Time of each voltage "
                      "step [default: %default]", default=1, type=float)
    parser.add_option("--repeats", dest="repeats",
                      help="how many times to run calculation for averaging"
                      " [default: %default]", default=1, type=int)
    parser.add_option("--dt", dest="dt", help="time step size"
                      " [default:%default]", default=0.1, type=float)
    parser.add_option("--file-name", dest="fileName", help="optional "
                      "output files name", default='')
    parser.add_option("--full", dest="fullOutput", help="if true the "
                      "results n and Q will be also saved [Default:%default]",
                      default=False, action='store_true')
    parser.add_option("-o", "--output-folder", dest="output_folder",
                      help="Output folder [default: current folder]",
                      default='.')

    # Disorder Parameters
    parser.add_option("--vg-avg", dest="VG_avg", help="Gate voltage average"
                      "[default: %default]", default=1, type=float)
    parser.add_option("--vg-std", dest="VG_std", help="Gate voltage std"
                      "[default: %default]", default=0, type=float)
    parser.add_option("--c-avg", dest="C_avg", help="capacitance of "
                      "junctions average [default: %default]",
                      default=1, type=float)
    parser.add_option("--c-std", dest="C_std", help="capacitance of "
                      "junctions std [default: %default]",
                      default=0, type=float)
    parser.add_option("--cg-avg", dest="CG_avg", help="Gate Capacitors "
                      "capacitance average [default: %default]",
                      default=1, type=float)
    parser.add_option("--cg-std", dest="CG_std", help="Gate Capacitors "
                      "capacitance std [default: %default]",
                      default=0, type=float)
    parser.add_option("--r-avg", dest="R_avg", help="junctions "
                      "resistance average [default: %default]",
                      default=1, type=float)
    parser.add_option("--r-std", dest="R_std", help="junctions "
                      "resistance std [default: %default]",
                      default=0, type=float)
    parser.add_option("--rg-avg", dest="RG_avg", help="Gate Resistors "
                      "resistance average [default: %default]",
                      default=1, type=float)
    parser.add_option("--rg-std", dest="RG_std", help="Gate Resistors "
                      "resistance std [default: %default]",
                      default=0, type=float)
    parser.add_option("--n-avg", dest="n0_avg", help="initial number of "
                      "electrons on each dot average [default:%default]",
                      default=0, type=float)
    parser.add_option("--n-std", dest="n0_std", help="initial number of "
                      "electrons on each dot std [default:%default]",
                      default=0, type=float)
    parser.add_option("--q-avg", dest="Q0_avg", help="initial charge on gate "
                      "capacitors average [default:%default]",
                      default=0, type=float)
    parser.add_option("--q-std", dest="Q0_std", help="initial charge on gate "
                      "capacitors std [default:%default]",
                      default=0, type=float)
    return parser.parse_args()

def saveParameters(path, fileName, options):

    optionsDict = options.__dict__
    with open(os.path.join(path,'runningParameters_' + fileName +
            ".txt"),mode='w') as f:
        for key in vars(options):
            f.write(key + " = " + str(optionsDict[key]) + "\n")

def create_random_array(M,N, avg, std, only_positive=False):
    if std == 0:
        return avg*np.ones((M, N))
    if only_positive:
        if avg == 0:
            shape = 1
            scale = 2
        else:
            shape = (avg/std)**2
            scale = std**2 / avg
        return 0.1 + np.random.gamma(shape=shape, scale=scale, size=(M,N))
    else:
        return np.random.normal(loc=avg, scale=std, size=(M,N))

if __name__ == "__main__":
    options, args = getOptions()
    columns = options.N
    Vext0 = options.Vmin
    VG = create_random_array(1, columns, options.VG_avg, options.VG_std,
                             False)
    dt = options.dt
    Q0 = create_random_array(1, columns, options.Q0_avg, options.Q0_std,
                             False)
    n0 = create_random_array(1, columns, options.n0_avg, options.n0_std, True)
    CG = create_random_array(1, columns, options.CG_avg, options.CG_std, True)
    RG = create_random_array(1, columns, options.RG_avg, options.RG_std, True)
    C = create_random_array(1, columns + 1, options.C_avg, options.C_std,
                            True)
    R = create_random_array(1, columns + 1, options.R_avg, options.R_std,
                            True)
    Vmax = options.Vmax
    Vstep = options.vStep
    tStep = options.tStep
    repeats = options.repeats
    savePath = options.output_folder
    fileName = options.fileName
    fullOutput = options.fullOutput

    if not os.path.exists(savePath):
        os.mkdir(savePath)
    elif not os.path.isdir(savePath):
        print("the given path exists but is a file")
        exit(0)

    saveParameters(savePath, fileName, options)
    run1DSimulation(Vext0, VG, dt, Q0, n0, CG, RG, C, R, columns,
                    Vmax, Vstep, tStep, repeats=repeats,
                    savePath=savePath, fileName=fileName, printState=False)