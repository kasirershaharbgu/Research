__author__ = 'shahar'

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from optparse import OptionParser
import os
from multiprocessing import Pool
from scipy.linalg import null_space
from scipy.integrate import cumtrapz
from scipy.signal import argrelextrema

MINIMUM_STEPS_PER_DOT = 1000

def flattenToColumn(a):
    return a.reshape((a.size, 1))

def flattenToRow(a):
    return a.reshape((1, a.size))

class NoElectronsOnDot(RuntimeError):
    pass

class DotArray:
    """
    An array of quantum dots connected to external voltage from left to
    right and to gate voltage
    """

    def __init__(self, rows, columns, VL, VR, VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, fastRelaxation=False):
        """
        Creates new array of quantum dots
        :param rows: number of rows (int)
        :param columns: number of columns (int)
        :param VL: left electrode voltage (double)
        :param VR: right electrode voltage (double)
        :param Q0: np array of initial charges (rowsXcolumns double array)
        :param n0: np array of initial electrons (rowsXcolumns double array)
        :param VG: voltage of the gates (rowsXcolumns double array)
        :param Ch: horizotal capacitances ((rowsXcolumns+1 double array)
        :param Cv: vertical capacitances (rows-1Xcolumns double array)
        :param Rh: vertical tunnelling ressistances (rowsXcolumns+1 double array)
        :param Rv: horizontal tunnelling ressistances (rowsXcolumns+1 double array)
        :param RG: ground resistances (rowsXcolumns double array)
        """
        self.fast_relaxation = fastRelaxation
        self.rows = rows
        self.columns = columns
        self.VL = VL
        self.VR = VR
        self.VG = VG
        self.Q = Q0
        self.n = n0
        self.CG = CG
        self.RG = RG
        self.Ch = Ch
        self.Cv = np.vstack((np.zeros((1,columns)),Cv)) if Cv.size else np.zeros((1,1))
        self.Rh = Rh
        self.Rv = Rv
        self.R = np.hstack((Rh.flatten(), Rh.flatten(), Rv.flatten(), Rv.flatten()))
        self.totalChargePassedRight = 0
        self.totalChargePassedLeft = 0
        self.no_tunneling_next_time = False
        self.createCapacitanceMatrix()
        self.setDiagonalizedJ()
        self.setConstWork()
        self.setConstMatrix()
        self.setConstNprimePart()
        self.saveCurrentMap = False
        self.Ih = None
        self.Iv = None

    def getRows(self):
        return self.rows

    def getColumns(self):
        return self.columns

    def resetCharge(self):
        """
        Reset the counter for total charge passed
        """
        self.totalChargePassedRight = 0
        self.totalChargePassedLeft = 0
        if self.saveCurrentMap:
            self.Ih = np.zeros(self.Rh.shape)
            self.Iv = np.zeros(self.Rv.shape)
        return True

    def changeVext(self, newVL, newVR):
        self.VL = newVL
        self.VR = newVR
        self.setConstWork()
        return True

    def getOccupation(self):
        return np.copy(self.n)

    def getGroundCharge(self):
        return np.copy(self.Q)

    def getCharge(self):
        return self.totalChargePassedRight, self.totalChargePassedLeft

    def getTimeStep(self):
        return self.timeStep

    def getCurrentMap(self):
        if self.saveCurrentMap:
            map = np.zeros((self.rows*2 - 1, self.columns+1))
            map[::2,:] = self.Ih
            map[1::2,:-1] = self.Iv
            return map
        else:
            raise NotImplementedError

    def currentMapOn(self):
        self.saveCurrentMap = True

    def currentMapOff(self):
        self.saveCurrentMap = False

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

    def setConstNprimePart(self):
        self._left_part_n_prime = np.copy(self.Ch[:, :-1])
        self._left_part_n_prime[:, 1:] = 0
        self._right_part_n_prime = np.copy(self.Ch[:, 1:])
        self._right_part_n_prime[:, :-1] = 0
        return True


    def getNprime(self):
        return flattenToColumn(self.n + self._left_part_n_prime*self.VL + self._right_part_n_prime*self.VR)

    def getbVector(self):
        """
        For developement of Q by dQ/dt=JQ +b
        """
        res = -self.invC.dot(self.getNprime()) + flattenToColumn(self.VG)
        return res / flattenToColumn(self.RG)

    def getJmatrix(self):
        """
        For developement of Q by dQ/dt=JQ +b
        """
        res = self.invC + np.diagflat(1/self.CG)
        return -res / np.repeat(flattenToColumn(self.RG),res.shape[1],axis=1)

    def setDiagonalizedJ(self):
        self._JeigenValues, self._JeigenVectors = np.linalg.eig(self.getJmatrix())
        self._JeigenValues = flattenToColumn(self._JeigenValues)
        self._JeigenVectorsInv = np.linalg.inv(self._JeigenVectors)
        self.timeStep = -10/np.min(self._JeigenValues)
        self.default_dt = -1/np.max(self._JeigenValues)
        if self.fast_relaxation:
            invMat = np.linalg.inv(self.invC + np.diagflat(1/self.CG))
            self._constQnPart = invMat.dot(flattenToColumn(self.VG))
            CGMat = np.repeat(1/flattenToRow(self.CG), invMat.shape[0], axis=0)
            self._matrixQnPart = invMat * CGMat - np.eye(invMat.shape[0])
        return True

    def developeQ(self, dt):
        Q0 = self._JeigenVectorsInv.dot(flattenToColumn(self.Q))
        b = self._JeigenVectorsInv.dot(self.getbVector())
        exponent = np.exp(self._JeigenValues*dt)
        self.Q = self._JeigenVectors.dot(Q0*exponent +(
             b/self._JeigenValues)*(exponent - 1)).reshape((self.rows, self.columns))
        return True

    def get_steady_Q_for_n(self):
        return self._constQnPart + self._matrixQnPart.dot(self.getNprime())

    def setConstWork(self):
        invCDiagMat = np.diag(np.copy(self.invC)).reshape((self.rows,self.columns))
        lowerCDiag = np.pad(np.diag(np.copy(self.invC),k=-1),((0,1),),mode='constant').reshape((self.rows,self.columns))
        lowerCDiag[:,-1] = 0
        lowerCDiag  = np.pad(lowerCDiag,((0,0),(1,0)),mode='constant')
        upperCDiag = np.pad(np.diag(np.copy(self.invC),k=1),((0,1),),mode='constant').reshape((self.rows,self.columns))
        upperCDiag[:,-1] = 0
        upperCDiag  = np.pad(upperCDiag,((0,0),(1,0)),mode='constant')
        commonHorz = np.pad(invCDiagMat,((0,0),(1,0)),mode='constant') + np.pad(invCDiagMat,((0,0),(0,1)),mode='constant')\
                        - lowerCDiag - upperCDiag

        lowerNCDiag = np.diag(np.copy(self.invC),k=-self.columns).reshape((self.rows-1,self.columns))
        upperNCDiag = np.diag(np.copy(self.invC),k=self.columns).reshape((self.rows-1,self.columns))
        commonVert = invCDiagMat[1:,:] + invCDiagMat[:-1,:] - lowerNCDiag - upperNCDiag

        additionalLeft = np.zeros((self.rows,self.columns+1))
        additionalLeft[:,0] = self.VL
        additionalLeft[:,-1] = -self.VR
        leftConstWork = (0.5*commonHorz + additionalLeft).flatten()

        additionalRight = -additionalLeft
        rightConstWork = (0.5*commonHorz + additionalRight).flatten()

        vertConstWork = (0.5*commonVert).flatten()

        self.constWork = np.hstack((rightConstWork, leftConstWork, vertConstWork, vertConstWork))
        return True

    def setConstMatrix(self):
        firstHorzMat = np.zeros(((self.columns + 1) * self.rows, self.columns * self.rows))
        secondHorzMat = np.zeros(firstHorzMat.shape)
        firstLocations = np.ones(((self.columns + 1) * self.rows,))
        secondLocations = np.ones(firstLocations.shape)
        firstLocations[self.columns::self.columns + 1] = 0
        secondLocations[0::self.columns + 1] = 0
        firstHorzMat[firstLocations.astype(np.bool), :] = np.copy(self.invC)
        secondHorzMat[secondLocations.astype(np.bool), :] = np.copy(self.invC)
        self.horizontalMatrix = firstHorzMat - secondHorzMat
        firstVertMat = np.copy(self.invC[self.columns:, :])
        secondVertMat = np.copy(self.invC[:-self.columns, :])
        self.verticalMatrix = firstVertMat - secondVertMat
        return True

    def getWork(self):
        q = self.getNprime() + flattenToColumn(self.Q)
        variableRightWork = (self.horizontalMatrix.dot(q)).flatten()
        variableDownWork = (self.verticalMatrix.dot(q)).flatten()
        variableWork = np.hstack((variableRightWork, -variableRightWork, variableDownWork, -variableDownWork))
        return variableWork + self.constWork

    def getRates(self):
        work = self.getWork()
        work[work > 0] = 0
        return -work / self.R

    def getTimeInterval(self, randomNumber):
        """
        Calculates dt using Gillespie algorithm
        :return:
        """
        if self.fast_relaxation:
            Q_copy = np.copy(self.Q)
            self.Q = self.get_steady_Q_for_n()
        rates = self.getRates()
        sum_rates = np.sum(rates)
        if sum_rates == 0:
            dt = self.default_dt
            self.no_tunneling_next_time = True
        else:
            dt = np.log(1/randomNumber)/sum_rates
        if self.fast_relaxation:
            self.Q = Q_copy
        return dt

    def nextStep(self, dt, randomNumber):
        self.developeQ(dt)
        if self.no_tunneling_next_time:
            self.no_tunneling_next_time = False
            return True
        rates = self.getRates()
        if (rates == 0).all():
            return True
        cumRates = np.cumsum(rates)
        actionInd = np.searchsorted(cumRates, randomNumber*cumRates[-1])
        self.executeAction(actionInd)
        return True

    def tunnel(self, fromDot, toDot):
        if fromDot[1] == self.columns:
            self.totalChargePassedRight -= 1
        elif fromDot[1] == -1:
            self.totalChargePassedLeft -= 1
        else:
            self.n[fromDot] -= 1
        if toDot[1] == self.columns:
            self.totalChargePassedRight += 1
        elif toDot[1] == -1:
            self.totalChargePassedLeft += 1
        else:
            self.n[toDot] += 1
        return True

    def printState(self):
        for i in range(self.rows):
            for j in range(self.columns):
                print("At dot (" + str(i) + ", " + str(j) + ") the state "
                      "is: n= " + str(self.n[i,j]) + " and QG = " + str(
                    self.Q[i,j]))

    def __str__(self):
        rows = "Rows: " + str(self.rows)
        columns = "Columns: " + str(self.columns)
        VG = "VG: " + str(self.VG)
        CG = "CG: " + str(self.CG)
        RG = "RG: " + str(self.RG)
        Ch = "Ch: " + str(self.Ch)
        Cv = "Cv: " + str(self.Cv)
        Rh = "Rh: " + str(self.Rh)
        Rv = "Rv: " + str(self.Rv)
        res = "----Array Parameters-----"
        for param in [rows, columns, VG, CG, RG, Ch, Cv, Rh, Rv]:
            res += "\n" + param
        return res

    def executeAction(self, ind):
        horzSize = self.rows*(self.columns+1)
        vertSize = (self.rows - 1)*self.columns
        if ind < horzSize: # tunnel right:
            fromDot = (ind//(self.columns+1), (ind%(self.columns+1))-1)
            toDot = (fromDot[0], fromDot[1] + 1)
            if self.saveCurrentMap:
                self.Ih[toDot] += 1
        elif ind < horzSize*2: # tunnel left
            ind -= horzSize
            fromDot = (ind//(self.columns+1), ind%(self.columns+1))
            toDot = (fromDot[0], fromDot[1] - 1)
            if self.saveCurrentMap:
                self.Ih[fromDot] -= 1
        elif ind < horzSize*2 + vertSize: # tunnel down
            ind -= horzSize*2
            fromDot = (ind//self.columns, ind%self.columns)
            toDot = (fromDot[0] + 1, fromDot[1])
            if self.saveCurrentMap:
                self.Iv[fromDot] += 1
        else: # tunnel up
            ind -= (horzSize*2 + vertSize)
            fromDot = ((ind//self.columns) + 1, ind%self.columns)
            toDot = (fromDot[0] - 1, fromDot[1])
            if self.saveCurrentMap:
                self.Iv[toDot] -= 1
        self.tunnel(fromDot, toDot)
        return fromDot, toDot

class Simulator:
    """
    Preforms Gillespie simulation on an array of quantum dots
    :param rows: number of rows in the array
    :param columns: number of columns in the array
    """
    def __init__(self, rows, columns, VL0, VR0, VG,
                 Q0, n0, CG, RG, Ch, Cv, Rh, Rv, fastRelaxation):
        self.dotArray = DotArray(rows, columns, VL0, VR0, VG, Q0, n0, CG, RG,
                                 Ch, Cv, Rh, Rv, fastRelaxation)
        self.VL = VL0
        self.VR = VR0
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

    def getArrayParameters(self):
        return str(self.dotArray)

    def executeStep(self, printState=False):
        r = np.random.rand(2)
        # for debug
        # r = next(self.randomGen)
        dt = self.dotArray.getTimeInterval(r[0])
        self.dotArray.nextStep(dt,r[1])
        if printState:
            self.printState()
        return dt

    def calcCurrent(self, t,print=False,fullOutput=False, currentMap=False):
        if currentMap:
            self.dotArray.currentMapOn()
        self.dotArray.resetCharge()
        curr_t = 0
        steps = 0
        if fullOutput:
            n_sum = np.zeros((self.dotArray.getRows(), self.dotArray.getColumns()))
            Q_sum = np.zeros((self.dotArray.getRows(), self.dotArray.getColumns()))
        min_steps = MINIMUM_STEPS_PER_DOT*(self.dotArray.getColumns()*self.dotArray.getRows())
        while curr_t < t or steps < min_steps:
            if fullOutput:
                curr_n = self.dotArray.getOccupation()
                curr_Q = self.dotArray.getGroundCharge()
            dt = self.executeStep(printState=print)
            steps += 1
            curr_t += dt
            if fullOutput:
                n_sum += curr_n*dt
                Q_sum += curr_Q*dt
        rightCurrent, leftCurrent = self.dotArray.getCharge()
        res = (rightCurrent/curr_t, -leftCurrent/curr_t)
        if fullOutput:
            res = res + (n_sum/curr_t, Q_sum/curr_t)
        if currentMap:
            res = res + (self.dotArray.getCurrentMap()/curr_t,)
            self.dotArray.currentMapOff()
        return res

    def checkSteadyState(self, rep):
        Ileft = []
        Iright = []
        t = self.dotArray.getTimeStep()
        for r in range(rep):
            # running once to get to steady state
            self.calcCurrent(t)
            # now we are in steady state calculate current
            rightCurrent, leftCurrnet = self.calcCurrent(t)
            Ileft.append(leftCurrnet)
            Iright.append(rightCurrent)
        x = t*np.arange(rep)
        print("Mean right: " + str(np.average(Iright))+ " Var right: " + str(np.var(Iright)))
        print("Mean left: " + str(np.average(Ileft))+ " Var left: " + str(np.var(Ileft)))
        print("Tail mean right: " + str(np.average(Iright[len(Iright)//4:])) + " Tail var right: " + str(np.var(Iright[len(Iright)//4:])))
        print("Tail mean left: " + str(np.average(Ileft[len(Ileft)//4:])) + " Tail var left: " + str(np.var(Ileft[len(Ileft)//4:])))
        plt.figure()
        plt.plot(x,Ileft, x,Iright, x, (np.array(Ileft) - np.array(Iright))**2)
        plt.figure()
        plt.plot((np.array(Ileft) - np.array(Iright)) ** 2)
        plt.show()

    def calcIV(self, Vmax, Vstep, fullOutput=False, print=False, currentMap=False):
        I = []
        if fullOutput:
            ns = []
            Qs = []
        if currentMap:
            Imaps = []
        tStep = self.dotArray.getTimeStep()
        VL_vec = np.arange(self.VL, Vmax+self.VR, Vstep)
        VL_vec = np.hstack((VL_vec, np.flip(VL_vec)))
        for VL in VL_vec:
            self.dotArray.changeVext(VL, self.VR)
            # running once to get to steady state
            self.calcCurrent(tStep)
            # now we are in steady state calculate current
            stepRes = self.calcCurrent(tStep, print=print, fullOutput=fullOutput, currentMap=currentMap)
            rightCurrent = stepRes[0]
            leftCurrnet = stepRes[1]
            if fullOutput:
                ns.append(stepRes[2])
                Qs.append(stepRes[3])
            if currentMap:
                Imaps.append(stepRes[-1])
            I.append((rightCurrent+leftCurrnet)/2)
        res = (np.array(I), VL_vec - self.VR)
        if fullOutput:
            res = res + (np.array(ns), np.array(Qs))
        if currentMap:
            res = res + (np.array(Imaps),)
        return res

    def printState(self):
        self.dotArray.printState()

class GraphSimulator:
    """
    Calculating staedy state current for an array of quantum dots by using linear programming
    """
    def __init__(self, rows, columns, VL0, VR0, VG,
                 Q0, n0, CG, RG, Ch, Cv, Rh, Rv, fastRelaxation):
        self.dotArray = DotArray(rows, columns, VL0, VR0, VG, Q0, n0, CG, RG,
                                 Ch, Cv, Rh, Rv, fastRelaxation)
        self.n0 = n0
        self.QG = Q0
        self.VL = VL0
        self.VR = VR0

        self.edgesMat = None
        self.states = None
        self.prob = None
        self.lyaponuv = None
        self.rates_diff_left = None
        self.rates_diff_right = None

        self.set_constant_matrices()
    def reshape_to_array(self,a):
        return a.reshape((self.dotArray.getRows(), self.dotArray.getColumns()))

    def getArrayParameters(self):
        return str(self.dotArray)


    def buildGraph(self, Q):
        states = []
        states_dict = dict()
        edges = []
        self.dotArray.Q = np.copy(Q)
        states.append(np.copy(self.n0).flatten())
        tup_n = tuple(self.n0.flatten())
        states_dict[tup_n] = 0
        left_rates_diff = []
        right_rates_diff = []
        current_state_ind = 0
        next_state_ind = 1
        while current_state_ind < len(states):
            self.dotArray.n = self.reshape_to_array(np.copy(states[current_state_ind]))
            edges_line = [0] * (current_state_ind + 1)
            rates = self.dotArray.getRates()
            for ind,rate in enumerate(rates):
                if rate > 0: # add new edge
                    fromDot, toDot = self.dotArray.executeAction(ind)
                    n = np.copy(self.dotArray.n).flatten()
                    tup_n = tuple(n)
                    if tup_n not in states_dict:
                        states_dict[tup_n] = next_state_ind
                        states.append(n)
                        next_state_ind += 1
                        edges_line.append(rate)
                    else:
                        edges_line[states_dict[tup_n]] += rate
                    self.dotArray.tunnel(toDot, fromDot)
            edges.append(edges_line)
            current_state_ind += 1
            left_diff, right_diff = self.get_edge_rates_diff(rates)
            left_rates_diff.append(left_diff)
            right_rates_diff.append(right_diff)
        edgesMat = np.zeros((len(edges),len(edges[-1])))
        for ind, line in enumerate(edges):
            edgesMat[ind,:len(line)] = line
        diagonal = np.sum(edgesMat,axis=1)
        self.edgesMat = edgesMat - np.diagflat(diagonal)
        self.states = np.array(states)
        self.rates_diff_left = np.array(right_rates_diff)
        self.rates_diff_right = np.array(left_rates_diff)

    def find_probabilities(self, Q):
        self.buildGraph(Q)
        # steady state probabilities are the belong to null-space of edgesMat.T
        null = null_space(self.edgesMat.T)
        sol = []
        for i in range(null.shape[1]):
            candidate = null[:,i]
            if (candidate >= 0).all() or (candidate <= 0).all():
                sol.append(candidate / np.sum(candidate))
        if len(sol) > 1:
            print("Warning - more than one solution for steady state probabilities")
        elif len(sol) < 1:
            self.prob = None
        else:
            self.prob = sol[0]

    def get_average_state(self, Q):
        self.find_probabilities(Q)
        average_state = np.sum(np.multiply(self.states.T, self.prob),axis=1)
        return self.reshape_to_array(average_state)

    def get_average_Qn(self, Q):
        self.dotArray.n = self.get_average_state(Q)
        return self._constQnPart + self._matrixQnPart.dot(self.dotArray.getNprime())

    def set_lyaponuv(self, Qmin, Qmax):
        Qmin = Qmin.flatten()
        Qmax = Qmax.flatten()
        coordinates = [np.arange(Qmin[i],Qmax[i],0.1) for i in range(Qmin.size)]
        grid = np.meshgrid(coordinates)
        grid_array = np.array(grid).T
        res = np.zeros(grid[0].shape)
        it = np.nditer(grid[0], flags=['multi_index'])
        while not it.finished:
            index = it.multi_index
            curr_Q = grid_array[index,:]
            diff_from_equi = curr_Q.flatten() - self.get_average_Qn(curr_Q.reshape((self.dotArray.getRows(),
                                     self.dotArray.getColumns())))
            res[it.multi_index] = diff_from_equi
            it.iternext()
        for axis in range(len(res.shape)):
            res = cumtrapz(res,grid[axis],axis=axis,initial=0)
        self.lyaponuv = res
        self.Q_grid = grid_array

    def find_next_QG(self):
        self.set_lyaponuv(self.QG - 1, self.QG + 1)
        peaks = argrelextrema(self.lyaponuv, np.less_equal)
        Qind = np.argmin(np.abs(self.Q_grid[peaks] - self.QG))
        self.QG = self.Q_grid[peaks][Qind]

    def calcCurrent(self, fullOutput=False):
        self.find_next_QG()
        self.dotArray.Q = self.reshape_to_array(self.QG)
        left_current = np.sum(self.prob*self.rates_diff_left.T)
        right_current = np.sum(self.prob*self.rates_diff_right.T)
        if fullOutput:
            n_avg = self.get_average_state(self.QG)
            return right_current, left_current, n_avg, QG
        return right_current, left_current

    def get_edge_rates_diff(self, rates):
        M = self.dotArray.rows
        N = self.dotArray.columns + 1
        horzSize = M * N
        right_tunneling_rates = rates[:horzSize].reshape((M, N))
        left_tunneling_rates = rates[horzSize:2*horzSize].reshape((M, N))
        left_diff = right_tunneling_rates[:,0] - left_tunneling_rates[:,0]
        right_diff = right_tunneling_rates[:,-1] - left_tunneling_rates[:,-1]
        return left_diff, right_diff

    def calcIV(self, Vmax, Vstep, fullOutput=False, print=False, currentMap=False):
        I = []
        if fullOutput:
            ns = []
            Qs = []
        VL_vec = np.arange(self.VL, Vmax + self.VR, Vstep)
        VL_vec = np.hstack((VL_vec, np.flip(VL_vec)))
        for VL in VL_vec:
            self.dotArray.changeVext(VL, self.VR)
            # running once to get to steady state
            self.calcCurrent()
            # now we are in steady state calculate current
            if fullOutput:
                rightCurrent, leftCurrnet, n, Q = self.calcCurrent(fullOutput=True)
                ns.append(n)
                Qs.append(Q)
            else:
                rightCurrent, leftCurrnet = self.calcCurrent()
            I.append((rightCurrent + leftCurrnet) / 2)
        if fullOutput:
            return np.array(I), VL_vec - self.VR, np.array(ns), np.array(Qs)
        return np.array(I), VL_vec - self.VR

def runSingleSimulation(VL0, VR0, VG0, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, rows, columns,
                        Vmax, Vstep, fullOutput=False, printState=False, useGraph=False, currentMap=False,
                        fast_relaxation=False):
    if useGraph:
        simulator = GraphSimulator(rows, columns, VL0, VR0, np.array(VG0),
                                    np.array(Q0), np.array(n0), np.array(CG),
                                    np.array(RG), np.array(Ch), np.array(Cv),
                                    np.array(Rh), np.array(Rv), fast_relaxation)
    else:
        simulator = Simulator(rows, columns, VL0, VR0, np.array(VG0),
                              np.array(Q0), np.array(n0), np.array(CG),
                              np.array(RG), np.array(Ch), np.array(Cv),
                              np.array(Rh), np.array(Rv), fast_relaxation)
    out = simulator.calcIV(Vmax, Vstep,
                           fullOutput=fullOutput, print=printState, currentMap=currentMap)
    array_params = simulator.getArrayParameters()
    return out + (array_params,)

def runFullSimulation(VL0, VR0, VG0, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, rows, columns,
                      Vmax, Vstep,repeats=1, savePath=".", fileName="", fullOutput=False,
                      printState=False, checkSteadyState=False, useGraph=False, fastRelaxation=False,
                      currentMap=False):
    if checkSteadyState:
        simulator = Simulator(rows, columns, VL0, VR0, np.array(VG0),
                              np.array(Q0), np.array(n0), np.array(CG),
                              np.array(RG), np.array(Ch), np.array(Cv),
                              np.array(Rh), np.array(Rv), fastRelaxation)
        simulator.checkSteadyState(repeats)
        exit(0)


    basePath = os.path.join(savePath, fileName)
    Is = []
    if fullOutput:
        ns = []
        Qs = []
    if currentMap:
        Imaps = []
    pool = Pool(processes=repeats)
    results = []
    for repeat in range(repeats):
        res = pool.apply_async(runSingleSimulation,
                                (VL0, VR0, VG0, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, rows, columns,
                                 Vmax, Vstep, fullOutput, printState, useGraph, currentMap, fastRelaxation))
        results.append(res)
    for res in results:
        result=res.get()
        I = result[0]
        currentMapInd = 2
        if fullOutput:
            n = result[2]
            Q = result[3]
            ns.append(n)
            Qs.append(Q)
            currentMapInd=4
        if currentMap:
            Imaps.append(result[currentMapInd])
        Is.append(I)
    V = result[1]
    params = result[-1]
    # dbg
    # results = []
    # for repeat in range(repeats):
    #     res = runSingleSimulation(VL0, VR0, VG0, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, rows, columns,
    #                              Vmax, Vstep, fullOutput, printState, useGraph, fast_relaxation)
    #     results.append(res)
    # for res in results:
    #     if fullOutput:
    #         I,V,n,Q,params = res
    #         ns.append(n)
    #         Qs.append(Q)
    #     else:
    #         I, V, params = res
    #     Is.append(I)
    # dbg

    avgI = np.mean(np.array(Is), axis=0)
    if fullOutput:
        avgN = np.mean(np.array(ns), axis=0)
        avgQ = np.mean(np.array(Qs), axis=0)
        np.save(basePath + "_n", avgN)
        np.save(basePath + "_Q", avgQ)
        np.save(basePath + "_full_I", np.array(Is))
    fig = plt.figure()
    plt.plot(V[:V.size//2], avgI[:I.size//2], '.b',V[V.size//2:],
             avgI[I.size//2:],
             '.r')
    plt.savefig(basePath + "_IV.png")
    np.save(basePath + "_I", avgI)
    np.save(basePath + "_V", V)
    plt.close(fig)
    if currentMap:
        avgImaps = np.mean(np.array(Imaps),axis=0)
        saveCurrentMaps(avgImaps,V, basePath + "_Imap")
    return params

def saveCurrentMaps(Imaps, V, path):
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=24, bitrate=1800)
    fig = plt.figure()
    Imax = np.max(Imaps)
    M,N = Imaps[0].shape
    im = plt.imshow(np.zeros(((M // 2) * 3 + 1, N * 3)), vmin=-Imax / 2,
                    vmax=Imax / 2, animated=True, cmap='plasma', aspect='equal')
    text = plt.text(1, 1, 'Vext = 0')
    plt.colorbar()
    frames = [(Imaps[i], V[i]) for i in range(len(V))]
    im_ani = animation.FuncAnimation(fig, plotCurrentMaps(im, text, M, N),
                                     frames=frames, interval=100,
                                     repeat_delay=1000,
                                     blit=True)
    im_ani.save(path + '.mp4', writer=writer)
    plt.close(fig)

def plotCurrentMaps(im, text, M, N):
    '''
    updating the plot to current currents map
    :return: image for animation
    '''
    J = np.zeros(((M//2)*3+1,N*3))
    vertRows = np.arange(0,(M//2)*3+1,3)
    vertCols = np.repeat(np.arange(0,3*N,3),2)
    vertCols[1::2] += 1
    horzRows = np.repeat(np.arange(1, (M//2)*3+1, 3),2)
    horzRows[1::2] += 1
    horzCols = np.arange(2,3*N,3)

    def updateCurrent(result):
        I, Vext = result
        if I is None:
            return im
        J[np.ix_(vertRows, vertCols)] = np.repeat(I[0:M:2,:],2,axis=1)
        J[np.ix_(horzRows, horzCols)] = np.repeat(I[1:M:2,:],2,axis=0)
        im.set_array(J)
        text.set_text('Vext = ' + str(Vext))
        return im,text
    return updateCurrent

def getOptions():
    parser = OptionParser(usage= "usage: %prog [options]")

    # Normal parameters
    parser.add_option("-M", "--height", dest="M", help="number of lines in "
                      "the array [default: %default]", default=1, type=int)
    parser.add_option("-N", "--width", dest="N", help="number of columns in "
                      "the array [default: %default]", default=1, type=int)
    parser.add_option("--vr", dest="VR", help="right electrode voltege [default: %default]", default=0, type=float)
    parser.add_option("--vmin", dest="Vmin", help="minimum external voltage"
                      " [default: %default]", default=0, type=float)
    parser.add_option("--vmax", dest="Vmax", help="maximum external voltage"
                      " [default: %default]", default=10, type=float)
    parser.add_option("--vstep", dest="vStep", help="size of voltage step["
                    "default: %default]", default=1, type=float)
    parser.add_option("--repeats", dest="repeats",
                      help="how many times to run calculation for averaging"
                      " [default: %default]", default=1, type=int)
    parser.add_option("--file-name", dest="fileName", help="optional "
                      "output files name", default='')
    parser.add_option("--distribution", dest="dist", help="probability distribution to use [Default:%default]",
                      default='uniform')
    parser.add_option("--full", dest="fullOutput", help="if true the "
                      "results n and Q will be also saved [Default:%default]",
                      default=False, action='store_true')
    parser.add_option("--graph", dest="use_graph", help="if true a simulation using graph solution for master equation"
                                                        "will be used [Default:%default]",
                      default=False, action='store_true')
    parser.add_option("--current-map", dest="current_map", help="if true  clip of current distribution during"
                                                                " simulation will be plotted [Default:%default]",
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

def saveParameters(path, fileName, options, array_params):
    optionsDict = options.__dict__
    with open(os.path.join(path,'runningParameters_' + fileName +
            ".txt"),mode='w') as f:
        f.write("-------Running parameters--------\n")
        for key in vars(options):
            f.write(key + " = " + str(optionsDict[key]) + "\n")
        f.write(array_params)

def create_random_array(M,N, avg, std, dist, only_positive=False):
    if std == 0:
        res = avg*np.ones((M, N))
    elif dist == 'uniform':
        res = np.random.uniform(low=avg-std, high=avg+std,size=(M,N))
    elif dist == 'two_points':
        if M == 1 and N == 2:
            return np.array([[avg-std,avg+std]])
        if M == 1 and N == 3:
            return np.array([[avg-std,avg+std, avg-std]])
        r = np.random.rand(M,N)
        r[r > 0.5] = avg + std
        r[r <= 0.5] = avg - std
        res = r
    elif dist == 'normal':
        res = np.random.normal(loc=avg, scale=std, size=(M, N))
    elif dist == 'gamma':
        if avg == 0:
            shape = 1
            scale = 2
        else:
            shape = (avg/std)**2
            scale = std**2 / avg
        res = 0.1 + np.random.gamma(shape=shape, scale=scale, size=(M,N))
    if only_positive and (res <= 0).any():
        print("Warnning, changing to positive distribution")
        res = res + 0.1 - np.min(res.flatten())
    return res


if __name__ == "__main__":
    # Initializing Running Parameters
    options, args = getOptions()
    rows = options.M
    columns = options.N
    VR0 = options.VR
    VL0 = VR0 + options.Vmin
    dist = options.dist
    VG = create_random_array(rows, columns, options.VG_avg, options.VG_std, dist,
                             False)
    Q0 = create_random_array(rows, columns, options.Q0_avg, options.Q0_std, dist,
                             False)
    n0 = create_random_array(rows, columns, options.n0_avg, options.n0_std, dist,False)
    CG = create_random_array(rows, columns, options.CG_avg, options.CG_std, dist,True)
    RG = create_random_array(rows, columns, options.RG_avg, options.RG_std, dist,True)
    Ch = create_random_array(rows, columns + 1, options.C_avg, options.C_std, dist,
                            True)
    Cv = create_random_array(rows - 1, columns, options.C_avg, options.C_std, dist,
                            True)
    Rh = create_random_array(rows, columns + 1, options.R_avg, options.R_std, dist,
                            True)
    Rv = create_random_array(rows - 1, columns, options.R_avg, options.R_std, dist,
                            True)
    Vmax = options.Vmax
    Vstep = options.vStep
    repeats = options.repeats
    savePath = options.output_folder
    fileName = options.fileName
    fullOutput = options.fullOutput
    use_graph = options.use_graph
    fast_relaxation = options.RG_avg * options.CG_avg < options.C_avg * options.R_avg
    current_map = options.current_map
    # Debug
    # rows = 1
    # columns = 1
    # VR0 = 0
    # VL0 = 0
    # VG = [[0]]
    # Q0 = [[0]]
    # n0 = [[0]]
    # CG = [[1]]
    # RG = [[1000]]
    # Ch = [[1,1]]
    # Cv = [[]]
    # Rh = [[1,10]]
    # Rv = [[]]
    # Vmax = 2
    # Vstep = 0.1
    # repeats = 1
    # savePath = "dbg_graph"
    # fileName = "dbg_graph"
    # fullOutput = False
    # use_graph = True

    # Running Simulation
    if not os.path.exists(savePath):
        os.mkdir(savePath)
    elif not os.path.isdir(savePath):
        print("the given path exists but is a file")
        exit(0)
    array_params = runFullSimulation(VL0, VR0, VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, rows,  columns,
                          Vmax, Vstep, repeats=repeats, savePath=savePath, fileName=fileName, fullOutput=fullOutput,
                          printState=False, useGraph=use_graph, fastRelaxation=fast_relaxation, currentMap=current_map)
    saveParameters(savePath, fileName, options, array_params)
    exit(0)
