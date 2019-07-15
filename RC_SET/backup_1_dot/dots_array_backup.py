__author__ = 'shahar'

import numpy as np
from matplotlib import pyplot as plt
from optparse import OptionParser
import os

class NoElectronsOnDot(RuntimeError):
    pass

class Dot:
    """
    A quantum dot connected to gate trough capacitor and resistor (RC-SET)
    """
    def __init__(self, Q0, n0, VG, CG, RG, VL, CL, RL, VR, CR, RR, VU=None,
                 CU=None, RU=None, VD=None, CD=None, RD=None):
        """
        Constructor for a single dot
        :param Q0: Initial charge on the gate capacitor (e units)
        :param n0: Initial number of electrons on the dot
        properties of connections:
        V - potential energy, C - capacitance, R - resistance
        the different possible connections are:
        G - gate, L- left, R- right, D - down, U - up
        """
        self.QG = Q0 # Charge on gate capacitor in units of e
        self.n = n0 # Number of electrons on the dot
        # connection properties:
        # Gate:
        self.VG = VG
        self.CG = CG
        self.RG = RG
        # Left:
        self.VL = VL
        self.CL = CL
        self.RL = RL
        # Right:
        self.VR = VR
        self.CR = CR
        self.RR = RR
        if VU:
            # Up: None if the dot is in the upper row in the array
            self.VU = VU
            self.CU = CU
            self.RU = RU
        if VD:
            # Down: None if the dot is in the upper row in the array
            self.VD = VD
            self.CD = CD
            self.RD = RD
        self.nChanged = True
        self._calcTau()


    def addElectron(self):
        """
        Adds electron to the dot
        """
        self.n += 1
        self.nChanged = True

    def removeElectron(self):
        """
        Removes electron from the dot
        """
        if self.n==0:
            raise NoElectronsOnDot
        self.n -= 1
        self.nChanged = True

    def changeVoltages(self, VL, VR, VU=None, VD=None):
        self.VL = VL
        self.VR = VR
        if VU:
            self.VU = VU
        if VD:
            self.VD = VD
        self._calcQn()


    def getNumberOfElectrons(self):
        """
        :return: the number of rlectrons on the dot
        """
        return self.n

    def getGateCharge(self):
        """
        :return: the charge on the gate capacitor
        """
        return self.QG

    def _calcTau(self):
        """
        Calculates relaxation time for the dot.
        """
        # TODO: upgrade to match dot that is connected from up and down
        self.tau = (self.RG * self.CG * (self.CL + self.CR)) / \
                   (self.CL + self.CR + self.CG)

    def _calcQn(self):
        """
        Calculating equilibrium value for the charge on the gate capacitor
        :return:
        """
        # TODO: upgrade to match dot that is connected from up and down
        self.Qn = self.tau * (((self.CL * self.VL +self.CR*self.VR + self.n) /
                  (self.RG*(self.CL + self.CR))) - self.VG/self.RG)


    def getRate(self, side, addElectron, includeVoltageWork=True):
        """
        Calculates the work done by the system when adding/removing
        electron trough the given side
        :param side: (R, L, U, D) the side trough the electron will be
        removed/added
        :param addElectron: If True an electron will be added, else an
        electron will be removed
        :return: The work done and the tunneling rate
        """
        # TODO: upgrade to match dot that is connected from up and down and
        #  split to energy difference and work by voltage difference
        if self.n < 1 and not addElectron:
            w = 0
            relevantR = 1
        else:
            sideSign = -1 if side== "L" else 1
            addSign = -1 if addElectron else 1
            relevantC = self.CR if side == "L" else self.CL
            relevantR = self.RL if side == "L" else self.RR
            if includeVoltageWork:
                w = (1/(self.CL + self.CR)) * \
                    (-0.5 + addSign * (2*(self.n - self.QG) +
                                    sideSign*relevantC*(self.VL - self.VR)))
            else:
                w = (1/(self.CL + self.CR)) * \
                    (-0.5 +  addSign *(self.n - self.QG))
        # TODO: upgrade for T > 0
        return w/relevantR

    def developQ(self, t):
        """
        Letting the system develop with no tunneling for time t
        :param t: time interval
        """
        if self.nChanged:
            self._calcQn()
            self.nChanged = False
        self.QG = (self.QG - self.Qn)*np.exp(-t/self.tau) + self.Qn


class DotArray:
    """
    An array of quantum dots connected to external voltage from left to
    right and to gate voltage
    """

    def __init__(self, rows, columns, Vext, VG, Q0, n0, CG, RG, CL,
                 RL, CR, RR, CU=None, RU=None, CD=None, RD=None):
        """
        Creates new array of quantum dots
        :param rows: number of rows
        :param columns: number of columns
        :param Vext: external voltage
        :param Q0: np array of initial charges
        :param n0: np array of initial electrons
        :param VG: voltage of the gate
        np arrays of other parameters of the dots:
        C - capacitance, R - resistance
        the different possible connections are:
        G - gate, L- left, R- right, D - down, U - up
        """
        self.rows = rows
        self.columns = columns
        self.Vext = Vext
        self.VG = VG
        VL, VR, _, _ = self._calcVoltages()
        self.array = [[None]*columns]*rows
        for i in range(rows):
            for j in range(columns):
                # TODO: upgrade for more than 1 row
                self.array[i][j] = Dot(Q0[i, j], n0[i, j], self.VG, CG[i, j],
                                       RG[i, j], VL[i, j], CL[i, j], RL[i, j],
                                       VR[i, j], CR[i, j], RR[i, j])
        self.totalChargePassed = 0
        # TODO: upgrade for more than 1 dot
        self.nearNeighbors = ((0,-1), (0,1))

    def _calcVoltages(self):
        """
        Calculating voltages for the array
        :return: np array of voltages (VL, VR, VU, VD)
        """
        # TODO: upgarde for more then 1 dot
        VL = np.array([[self.Vext/2]])
        VR = np.array([[-self.Vext/2]])
        return VL, VR, None, None

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
        self.Vext = newVext
        VL, VR, _, _ = self._calcVoltages()
        for i in range(self.rows):
            for j in range(self.columns):
                # TODO: upgrade for more than 1 row
                self.array[i][j].changeVoltages(VL[i][j], VR[i][j])

    def getCharge(self):
        return self.totalChargePassed

    def getRate(self, fromDot, toDot):
        """
        Returns the work done when tunneling from dot
        :param fromDot: index of the dot from where the electron tunnels ,
        column -1 means left electrode
        :param toDot: index of the dot to where to tunnel, column index =
        number of columns means right electrode
        :return: the work done by this tunneling event
        """
        # TODO: upgrade for more than 1 dot
        if fromDot[1] < 0:
            return self.array[toDot[0]][toDot[1]].getRate("L", True)
        elif fromDot[1] >= self.rows:
            return self.array[toDot[0]][toDot[1]].getRate("R", True)
        elif toDot[1] < 0:
            return self.array[fromDot[0]][fromDot[1]].getRate("L", False)
        elif toDot[1] >= self.rows:
            return self.array[fromDot[0]][fromDot[1]].getRate("R", False)
        else:
            raise(NotImplementedError)


    def getAllPossibleRates(self):
        """
        :return: dictionary with keys-> possible tunnels and items work done
        by this tunnel
        """
        # TODO: upgrade for more than 1 dot
        rates = []
        for row in range(0, self.rows,2,):
            for col in range(0, self.columns, 2):
                rate = []
                for p in self.nearNeighbors:
                    current = (row, col)
                    neighbor = (row + p[0],col + p[1])
                    rate.append(self.getRate(current, neighbor))
                    rate.append(self.getRate(neighbor, current))
                rates.append(rate)
        return np.array(rates)

    def develope(self, t):
        for row in range(self.rows):
            for column in range(self.columns):
                self.array[row][column].developQ(t)

    def tunnel(self, fromDot, toDot):
        if self.columns > fromDot[1] >= 0:
            self.array[fromDot[0]][fromDot[1]].removeElectron()
        if self.columns > toDot[1] >= 0:
            self.array[toDot[0]][toDot[1]].addElectron()
        if toDot[1] == self.columns:
            self.totalChargePassed += 1
        elif fromDot[1] == self.columns:
            self.totalChargePassed -= 1

    def printState(self):
        for i in range(self.rows):
            for j in range(self.columns):
                n = self.array[i][j].getNumberOfElectrons()
                Q = self.array[i][j].getGateCharge()
                print("At dot (" + str(i) + ", " + str(j) + ") the state "
                      "is: n= " + str(n) + " and QG = " + str(Q))

    def getIndices(self, ind):
            return ind // self.columns, ind % self.columns

MINPROBLOG = np.log(1/0.7)
class Simulator:
    """
    Preforms Markov-Chain Monte-Carlo simulation on an array of quantum dots
    :param rows: number of rows in the array
    :param columns: number of columns in the array
    """
    def __init__(self, rows, columns, Vext0, VG0, dt,
                 Q0, n0, CG, RG, CL, RL, CR, RR):
        # TODO upgrade for more than 1 row
        self.dotArray = DotArray(rows, columns, Vext0, VG0, Q0, n0, CG, RG,
                                 CL, RL, CR, RR)
        self.t = 0
        self.dt = dt
        self.prob = None
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


    def _calcProbabilities(self):
        rates = self.dotArray.getAllPossibleRates()
        self.prob = np.zeros(rates.shape)
        relevant = rates > 0
        self.prob[relevant] = 1 - np.exp(-rates[relevant]* self.dt)
        self.prob = np.cumsum(self.prob,axis = 1)

    def executeStep(self, printState=False):
        self.dotArray.develope(self.dt)
        self._calcProbabilities()

        r = np.random.ranf((self.prob.shape[0],)) # TODO return
        # r = next(self.randomGen)

        tunnelToLeft, = np.where(r < self.prob[:,0])
        tunnelFromLeft, = np.where((r < self.prob[:,1]) ^ (r < self.prob[:,0]))
        tunnelToRight, = np.where((r < self.prob[:,2]) ^ (r < self.prob[:,1]))
        tunnelFromRight, = np.where((r < self.prob[:,3]) ^ (r < self.prob[:,2]))
        try:
            for tunnel in tunnelToLeft:
                row, col = self.dotArray.getIndices(tunnel)
                self.dotArray.tunnel((row, col), (row, col-1))
            for tunnel in tunnelToRight:
                row, col = self.dotArray.getIndices(tunnel)
                self.dotArray.tunnel((row, col), (row, col+1))
            for tunnel in tunnelFromLeft:
                row, col = self.dotArray.getIndices(tunnel)
                self.dotArray.tunnel((row, col - 1), (row, col))
            for tunnel in tunnelFromRight:
                row, col = self.dotArray.getIndices(tunnel)
                self.dotArray.tunnel((row, col + 1), (row, col))
            self.t += self.dt
            if printState:
                self.printState()
        except NoElectronsOnDot as e:
            self.dotArray.printState()
            print(self.prob)
            raise e

    def calcCurrent(self, t, fullOutput=False):
        # TODO: for now t must be multiplication of dt by whole number
        if fullOutput:
            n = []
            Q = []
        for steps in range(int(t // self.dt)):
            self.executeStep(printState=False)

            if fullOutput:
                n.append(self.dotArray.array[0][0].getNumberOfElectrons())
                Q.append(self.dotArray.array[0][0].getGateCharge())

        current = self.dotArray.getCharge() / t
        if fullOutput:
            return current,n, Q
        return current

    def calcIV(self, Vmax, Vstep, tStep, fullOutput=False):
        V = []
        I = []
        if fullOutput:
            Q = []
            n = []
        V0 = self.Vext
        while(self.Vext < Vmax):
            if fullOutput:
                current, stepn, stepQ = self.calcCurrent(tStep,
                                                         fullOutput=True)
                I.append(current)
                n.extend(stepn)
                Q.extend(stepQ)
            else:
                I.append(self.calcCurrent(tStep, fullOutput=False))
            V.append(self.Vext)
            self.Vext += Vstep
            self.dotArray.changeVext(self.Vext)
            self.dotArray.resetCharge()
        while(self.Vext > V0):
            if fullOutput:
                current, stepn, stepQ = self.calcCurrent(tStep,
                                                         fullOutput=True)
                I.append(current)
                n.extend(stepn)
                Q.extend(stepQ)
            else:
                I.append(self.calcCurrent(tStep, fullOutput=False))
            V.append(self.Vext)
            if count < 10:
                count +=1
            else:
                self.Vext -= Vstep
                self.dotArray.changeVext(self.Vext)
                self.dotArray.resetCharge()
        if fullOutput:
            return np.array(I), np.array(V), np.array(n), np.array(Q)
        return np.array(I), np.array(V)

    def printState(self):
        print("At t = " + str(self.t) + ":")
        self.dotArray.printState()


def runSingleDotSimulation(Vext0, VG0, dt, Q0, n0, CG, RG, CL, RL, CR, RR,
                           Vmax, Vstep, tStep,repeats=1, savePath=".",
                           fileName="", fullOutput=False):
    basePath = os.path.join(savePath, fileName)
    Is = []

    for repeat in range(repeats):
        simulator = Simulator(1, 1, Vext0, VG0, dt, np.array([[Q0]]),
                              np.array([[n0]]), np.array([[CG]]),
                              np.array([[RG]]), np.array([[CL]]),
                              np.array([[RL]]), np.array([[CR]]),
                              np.array([[RR]]))
        out = simulator.calcIV(Vmax, Vstep, tStep,
                                fullOutput=fullOutput)
        I = out[0]
        V = out[1]
        if fullOutput:
            n = out[2]
            Q = out[3]
        Is.append(I)
    avgI = np.mean(np.array(Is), axis=0)
    if fullOutput:
        fig = plt.figure()
        plt.plot(n, '.')
        plt.savefig(basePath + "_n.png")
        np.save(basePath + "_n.bin", n)
        plt.close(fig)

        fig = plt.figure()
        plt.plot(Q, '.')
        plt.savefig(basePath + "_Q.png")
        np.save(basePath + "_Q.bin", Q)
        plt.close(fig)

    fig = plt.figure()
    plt.plot(V[:V.size//2], avgI[:I.size//2], '.b',V[V.size//2:],
             avgI[I.size//2:],
             '.r')

    plt.savefig(basePath + "_IV.png")
    np.save(basePath + "_I.bin", avgI)
    np.save(basePath + "_V.bin", V)
    plt.close(fig)
    if fullOutput:
        return I, V, n, Q
    else:
        return I, V

def getOptions():
    parser = OptionParser(usage= "usage: %prog [options]")

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
    parser.add_option("--gvoltage", dest="VG", help="Gate voltage "
                      "[default: %default]", default=0, type=float)
    parser.add_option("--rightc", dest="CR", help="Right junction "
                      "capacitance[default: %default]",
                      default=1, type=float)
    parser.add_option("--leftc", dest="CL", help="Left junction "
                      "capacitance[default: %default]",
                      default=1, type=float)
    parser.add_option("--groundc", dest="CG", help="Gate Capacitor "
                      "capacitance[default: %default]",
                      default=1, type=float)
    parser.add_option("--rightr", dest="RR", help="Right junction "
                      "resistance[default: %default]",
                      default=1, type=float)
    parser.add_option("--leftr", dest="RL", help="Left junction "
                      "resistance [default: %default]",
                      default=1, type=float)
    parser.add_option("--groundr", dest="RG", help="Gate Resistor "
                      "resistance[default: %default]",
                      default=1, type=float)
    parser.add_option("--repeats", dest="repeats",
                      help="how many times to run calculation for averaging"
                      " [default: %default]", default=1, type=int)
    parser.add_option("-n", dest="n0", help="initial number of electrons on"
                      " dot [default:%default]", default=0, type=float)
    parser.add_option("-Q", dest="Q0", help="initial charge on gate capacitor"
                      " [default:%default]", default=0, type=float)
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

    return parser.parse_args()

def saveParameters(path, fileName, options):

    optionsDict = options.__dict__
    with open(os.path.join(path,'runningParameters_' + fileName +
            ".txt"),mode='w') as f:
        for key in vars(options):
            f.write(key + " = " + str(optionsDict[key]) + "\n")




if __name__ == "__main__":
    # options, args = getOptions()

    # Vext0 = options.Vmin
    # VG = options.VG
    # dt = options.dt
    # Q0 = options.Q0
    # n0 = options.n0
    # CG = options.CG
    # RG = options.RG
    # CL = options.CL
    # RL = options.RL
    # CR = options.CR
    # RR = options.RR
    # Vmax = options.Vmax
    # Vstep = options.vStep
    # tStep = options.tStep
    # repeats = options.repeats
    # savePath = options.output_folder
    # fileName = options.fileName
    # fullOutput = options.fullOutput

    # saveParameters(savePath, fileName, options)

    Q0 = 0
    n0 = 0
    CL = 1
    RL = 1
    Vmax = 5
    Vstep = 0.1
    tStep = 100
    repeats = 10
    fullOutput = False

    Vext0 = 2
    dt = 0.01

    VGs = [0]
    CGs = [1]
    RGs = [1000]
    CRs = [1]
    RRs = [3]

    for VG in VGs:
        for CG in CGs:
            for RG in RGs:
                for CR in CRs:
                    for RR in RRs:
                        savePath = "dbg" + str(
                            VG) + \
                                   "_CG_" + str(CG) + "_RG_" + \
                                   str(RG)
                        fileName = "CR_" + str(CR) + "_RR_"+ str(RR)
                        if not os.path.exists(savePath):
                            os.mkdir(savePath)


                        runSingleDotSimulation(Vext0, VG, dt, Q0, n0, CG, RG, CL, RL, CR, RR,
                           Vmax, Vstep, tStep, repeats=repeats,
                           savePath=savePath, fileName=fileName,
                           fullOutput=fullOutput)