__author__ = 'shahar'

import numpy as np
from matplotlib import pyplot as plt
from optparse import OptionParser
import os
from ast import literal_eval

def sumToColumnVector(a, axis):
    return a.sum(axis).reshape((a.shape[1+axis],1))

def padAndShift(a):
    tempZeros = np.zeros((a.shape[0],a.shape[1] + 2))
    aShiftedRight = tempZeros.copy()
    aShiftedLeft = tempZeros.copy()
    aPadded = tempZeros.copy()
    aShiftedRight[:,2:] = a
    aShiftedLeft[:,:-2] = a
    aPadded[:,1:-1] = a
    return aShiftedRight, aShiftedLeft, aPadded

class NoElectronsOnDot(RuntimeError):
    pass

class DotArray:
    """
    An array of quantum dots connected to external voltage from left to
    right and to gate voltage
    """

    def __init__(self, rows, columns, Vext, VG, Q0, n0, CG, RG, C,
                 R):
        """
        Creates new array of quantum dots
        :param rows: number of rows
        :param columns: number of columns
        :param Vext: external voltage
        :param Q0: np array of initial charges
        :param n0: np array of initial electrons
        :param VG: voltage of the gates (NXM array)
        np arrays of other parameters of the dots:
        C - capacitance, R - resistance
        the different possible connections are:
        G - gate, L- left, R- right, D - down, U - up
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
        self.C = C
        self.R = R
        self.totalChargePassed = 0
        self.createConstantMatrices()

    def getRows(self):
        return self.rows

    def getColumns(self):
        return self.columns

    def resetCharge(self):
        """
        Reset the counter for total charge passed
        """
        self.totalChargePassed = 0

    def changeVext(self, newVext):
        self.VL = newVext/2
        self.VR = -newVext/2
        self.createConstantMatrices()

    def getCharge(self):
        return self.totalChargePassed

    def createConstantMatrices(self):
        invC = 1/self.C
        invCMatrix = np.repeat(invC[:,:-1],self.columns,0)
        A = np.tril(invCMatrix)
        self._invCsqrt = np.sqrt(invC)
        self._invCG = 1/self.CG
        self._invRG = 1/self.RG
        self._invR = 1/self.R
        B = np.triu(np.ones((self.columns, self.columns)))
        M = np.linalg.inv(np.eye(self.columns) + self.C[0,-1]*invCMatrix)
        AdotM = A.dot(M)
        constant = self.C[0,-1] * (self.VR - self.VL)
        self._AdotMdotB = AdotM.dot(B)
        J = -np.diag(self._invRG.flatten()).dot(self._AdotMdotB + np.diag(
            self._invCG.flatten()))
        self._JeigenValues, self._JeigenVectors = np.linalg.eig(J)
        self._JeigenValues = self._JeigenValues.reshape((self._JeigenValues.shape[0], 1))
        self._JeigenVectorsInv = np.linalg.inv(self._JeigenVectors)
        self._bBase = -self._invRG.T*(constant *
                             sumToColumnVector(AdotM, -1) - self.VG.T + self.VL)
        Mtilde = np.diag(self._invCsqrt[0,:-1]).dot(M)
        alpha = Mtilde.dot(B)
        alphaAddition = -np.sqrt(self.C[0,-1])*A[-1,:].dot(alpha)
        self._alpha = np.vstack((alpha, alphaAddition))
        beta = constant*sumToColumnVector(Mtilde,-1)
        betaAddition = np.sqrt(self.C[0,-1])*(self.VR-self.VL - A[-1,
                                                                :].dot(beta))
        self._beta = np.vstack((beta, betaAddition))
        alphaShiftedRight, alphaShiftedLeft, alphaPadded = padAndShift(
            self._alpha)

        dalphaRight = alphaShiftedLeft - alphaPadded
        dalphaLeft = alphaShiftedRight - alphaPadded
        dalphaRightsqr = sumToColumnVector(dalphaRight*dalphaRight, 0)
        dalphaLeftsqr = sumToColumnVector(dalphaLeft*dalphaLeft, 0)

        self._deltaERightBase = 0.5*(dalphaRightsqr + dalphaRight.T.dot(
            self._beta))

        self._deltaELeftBase = 0.5*(dalphaLeftsqr + dalphaLeft.T.dot(
            self._beta))

    def getbVector(self):
        return self._JeigenVectorsInv.dot(self._bBase - (
            self._invRG.T*self._AdotMdotB.dot(
            self.n.T)))

    def developeQ(self, dt):
            Q0 = self._JeigenVectorsInv.dot(self.Q.T)
            b = self.getbVector()
            exponent = np.exp(self._JeigenValues*dt)
            self.Q = self._JeigenVectors.dot(Q0*exponent +(
                b/self._JeigenValues)*(exponent - 1)).T

    def getprobabilities(self, dt):
        alphaShiftedRight, alphaShiftedLeft, alphaPadded = padAndShift(
            self._alpha)
        dalphaRight = alphaShiftedLeft - alphaPadded
        dalphaLeft = alphaShiftedRight - alphaPadded


        commonElement = self._alpha.dot((self.Q+self.n).T)
        deltaERight = self._deltaERightBase + dalphaRight.T.dot(commonElement)
        deltaELeft = self._deltaELeftBase + dalphaLeft.T.dot(commonElement)


        Qjunctions = commonElement + self._beta
        QCratio = Qjunctions*self._invCsqrt.T
        tempZeros = np.zeros((QCratio.shape[0]+1,QCratio.shape[1]))
        deltaVRight = tempZeros.copy()
        deltaVLeft = tempZeros.copy()

        deltaVRight[:-1,:] = QCratio
        deltaVLeft[1:,:] = -QCratio

        tempZeros = np.zeros((self._invR.shape[0],self._invR.shape[1]+1))
        invRShiftedRight = tempZeros.copy()
        invRShiftedLeft = tempZeros.copy()

        invRShiftedRight[:,1:] = self._invR
        invRShiftedLeft[:,:-1] = self._invR

        rateRight = (deltaERight + deltaVRight).T* invRShiftedLeft
        rateLeft = (deltaELeft + deltaVLeft).T* invRShiftedRight
        noElectrons = np.hstack(([[False]], self.n == 0, [[False]]))
        rateRight[noElectrons] = 0
        rateLeft[noElectrons] = 0


        ratesCombined = np.vstack((rateRight, rateLeft)).flatten()
        prob = np.zeros(ratesCombined.shape)
        prob[ratesCombined < 0] = 1 - np.exp(ratesCombined[ratesCombined <
                                                         0]*dt)
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
        if (self.n < 0).any():
            raise NoElectronsOnDot

    def printState(self):
        for i in range(self.rows):
            for j in range(self.columns):
                print("At dot (" + str(i) + ", " + str(j) + ") the state "
                      "is: n= " + str(self.n[i,j]) + " and QG = " + str(
                    self.Q[i,j]))

    def executeStep(self, ind):
            side = ind // ((self.columns+2)*self.rows)
            dotIndex = ind % ((self.columns+2)*self.rows)
            dotRow = dotIndex // (self.columns+2)
            dotColumn = dotIndex % (self.columns+2) -1
            fromDot = (dotRow, dotColumn)
            if side == 0:
                toDot = (dotRow, dotColumn+1)
            if side == 1:
                toDot = (dotRow, dotColumn-1)
            self.tunnel(fromDot, toDot)

MINPROBLOG = np.log(1/0.7)
class Simulator:
    """
    Preforms Markov-Chain Monte-Carlo simulation on an array of quantum dots
    :param rows: number of rows in the array
    :param columns: number of columns in the array
    """
    def __init__(self, rows, columns, Vext0, VG0, dt,
                 Q0, n0, CG, RG, C, R):
        # TODO upgrade for more than 1 row
        self.dotArray = DotArray(rows, columns, Vext0, VG0, Q0, n0, CG, RG,
                                 C, R)
        self.t = 0
        self.dt = dt
        self.Vext = Vext0
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
        # TODO: for now t must be multiplication of dt by whole number
        if fullOutput:
            n = []
            Q = []
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


def run1DSimulation(Vext0, VG0, dt, Q0, n0, CG, RG, C, R, columns,
                    Vmax, Vstep, tStep,repeats=1, savePath=".",
                    fileName="", fullOutput=False, printState=False):
    if not os.path.exists(savePath):
        os.mkdir(savePath)
    elif not os.path.isdir(savePath):
        print("the given path exists but is a file")
        exit(0)
    basePath = os.path.join(savePath, fileName)
    Is = []

    for repeat in range(repeats):
        simulator = Simulator(1, columns, Vext0, np.array(VG0), dt,
                              np.array(Q0),
                              np.array(n0), np.array(CG),
                              np.array(RG), np.array(C),
                              np.array(R))
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
        return np.random.gamma(shape=shape, scale=scale, size=(M,N))
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

    saveParameters(savePath, fileName, options)
    run1DSimulation(Vext0, VG, dt, Q0, n0, CG, RG, C, R, columns,
                    Vmax, Vstep, tStep, repeats=repeats,
                    savePath=savePath, fileName=fileName, printState=False)