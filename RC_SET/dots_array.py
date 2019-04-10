__author__ = 'shahar'

import numpy as np


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
    def addElectron(self):
        """
        Adds electron to the dot
        """
        self.n += 1

    def removeElectron(self):
        """
        Removes electron from the dot
        """
        if self.n==0:
            raise("no electrons in dot")
        self.n -= 1

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
        self.Qn = self.tau * ((self.CL * (self.VL - self.VR) + self.n) /
                  (self.RG*(self.CL + self.CR)) - self.VG/self.RG)


    def getWork(self, side, addElectron):
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
        if self.n ==0 and not addElectron:
            w = np.inf
            relevantR = 1
        else:
            sideSign = -1 if side== "L" else 1
            addSign = -1 if addElectron else 1
            relevantC = self.CR if side == "L" else self.CL
            relevantR = self.RL if side == "L" else self.RR
            w = (1/(self.CL + self.CR)) * \
                (0.5 + addSign*(self.n - self.QG +
                                sideSign*relevantC*(self.VL - self.VR)))
        # TODO: upgrade for T > 0
        return w, w/relevantR

    def developeQ(self, t):
        """
        Letting the system develop with no tunneling for time t
        :param t: time interval
        """
        if self.nChanged:
            self._calcQn()
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
        self.possibleTunnels = (((0,-1),(0,0)), ((0,1), (0,0)),
                                ((0,0), (0, -1)), ((0,0), (0,1)))

    def _calcVoltages(self):
        """
        Calculating voltages for the array
        :return: np array of voltages (VL, VR, VU, VD)
        """
        # TODO: upgarde for more then 1 dot
        VL = np.array([[self.Vext]])
        VR = np.array([[0]])
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

    def getCharge(self):
        return self.totalChargePassed

    def getWork(self, fromDot, toDot):
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
            return self.array[toDot[0]][toDot[1]].getWork("L", True)
        elif fromDot[1] > self.rows:
            return self.array[toDot[0]][toDot[1]].getWork("R", True)
        elif toDot[1] < 0:
            return self.array[fromDot[0]][fromDot[1]].getWork("L", False)
        elif toDot[1] > self.rows:
            return self.array[fromDot[0]][fromDot[1]].getWork("R", False)
        else:
            # Not implemented yet
            return None

    def getAllPossibleWorks(self):
        """
        :return: dictionary with keys-> possible tunnels and items work done
        by this tunnel
        """
        # TODO: upgrade for more than 1 dot
        works = []
        rates = []
        for p in self.possibleTunnels:
            work, rate = self.getWork(p[0], p[1])
            works.append(work)
            rates.append(rate)
        return np.array(works), np.array(rates)

    def develope(self, t):
        for row in range(self.rows):
            for column in range(self.columns):
                self.array[row][column].developQ(t)

    def tunnel(self, index):
        event = self.possibleTunnels[index]
        fromDot = event[0]
        toDot = event[1]
        self.array[fromDot[0]][fromDot[1]].removeElectron()
        self.array[toDot[0]][toDot[1]].addElectron()

MINPROBLOG = np.log(1/0.7)
class Simulator:
    """
    Preforms Markov-Chain Monte-Carlo simulation on an array of quantum dots
    :param rows: number of rows in the array
    :param columns: number of columns in the array
    """
    def __init__(self, rows, columns, Vext0, VG0, dt):
        # TODO upgrade for more than 1 row
        Q0, n0, CG, RG, CL, RL, CR, RR = self._getParameters()
        self.dotArray = DotArray(rows, columns, Vext0, VG0, Q0, n0, CG, RG,
                                 CL,
                 RL, CR, RR)
        self.t = 0
        self.dt = dt
        self.prob = None

    def _getParameters(self):
        # TODO upgrade
        zero = np.zeros((1,1))
        one = np.ones((1,1))
        return zero, zero, one, one, one, one, one, one


    def _calcProbabilities(self):
        works, rates = self.dotArray.getAllPossibleWorks()
        self.prob = np.zeros(works.shape)
        # maxInd = np.argmax(works)
        # dt = MINPROBLOG / (works[maxInd] * rates[maxInd])
        relevant = works > 0
        self.prob[relevant] = 1 - np.exp(-works[relevant] * rates[relevant]
                                         * self.dt)

    def executeStep(self):
        self.dotArray.develope(self.dt)
        self._calcProbabilities()
        r = np.random.ranf(self.prob.shape)
        execute, = np.where(r < self.prob)
        # TODO change, good for one dot
        if execute.size > 1:
            execute = [np.argmax(self.prob)]
        for tunnel in execute:
            self.dotArray.tunnel(tunnel)
        self.t += self.dt

    def calcCurrent(self, t):
        # TODO: for now t must be multiplication of dt by whole number
        for steps in range(t // self.dt):
           self.executeStep()

        current = self.dotArray.getCharge() / t
        return current













