__author__ = 'shahar'

import numpy as np
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
from matplotlib import pyplot as plt
import matplotlib.animation as animation
# from mayavi import mlab
from optparse import OptionParser
import os
from multiprocessing import Pool
from scipy.linalg import null_space
from scipy.integrate import cumtrapz, dblquad
from scipy.interpolate import interp1d
from scipy.stats import norm
from mpl_toolkits.mplot3d import Axes3D
from ast import literal_eval
from copy import copy

# Gillespie Constants
MINIMUM_STEPS_PER_DOT = 1000
# Lyaponuv Constants
EPS = 0.0001
DQ = 0.1
Q_SHIFT = 2
GRAD_REP = 5
INI_LR = 0.001



def flattenToColumn(a):
    return a.reshape((a.size, 1))

def flattenToRow(a):
    return a.reshape((1, a.size))

def local_min_func(input):
    mid = len(input)//2
    return (input[mid] <= input).all()

def detect_local_minima(arr):
    # https://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the local minimum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood minimum, 0 otherwise)
    """
    neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
    return np.where(filters.generic_filter(arr, local_min_func, footprint=neighborhood,
                                           mode='constant', cval=np.min(arr)-1))
def simple_gadient_descent(grad, x0, eps=1e-10, lr=1e-3, max_iter=1e6, plot_lc=False):
    x=x0.flatten()
    curr_grad = grad(x)
    if plot_lc:
        gradvec = []
    iter=0
    while np.max(np.abs(curr_grad)) > eps and iter < max_iter:
        if plot_lc:
            gradvec.append(curr_grad)
        x = x - curr_grad*lr
        curr_grad = grad(x)
        iter += 1
    if plot_lc:
        plt.figure()
        plt.plot(gradvec)
    return x, iter < max_iter

# Methods for calculating superconducting related tunneling rates

def high_impadance_p(Ec,T,kappa):
    sigma = 4*(kappa**2)*Ec*T
    mu = kappa**2*Ec
    def f(x):
        return norm.pdf(x, loc=mu, scale=np.sqrt(sigma))
    return f

def fermi_dirac_dist(T):
    def f(x):
        return 1/(1 + np.exp(x/T))
    return f

def qp_density_of_states(energy_gap):
    def f(x):
        return np.abs(x) / np.sqrt(x**2-energy_gap**2)

    return f

def cp_tunneling(Ej, Ec, T):
    p = high_impadance_p(Ec, T, 2)
    def f(x):
        return (np.pi*Ej)**2 * p(x)
    return f

def qp_tunneling(Ec, gap, T):
    p = high_impadance_p(Ec, T, 1)
    N = qp_density_of_states(gap - EPS)
    fer = fermi_dirac_dist(T)

    def integrand(x2, x1, deltaE):
        return N(x1)*N(x2)*fer(-x1)*fer(x2)*p(x2-x1 + deltaE)

    def integrate_single(deltaE):
        part1 = dblquad(integrand, -np.infty, -gap, lambda x: -np.infty,
                                lambda x: -gap, args=(deltaE,))[0]
        part2 = dblquad(integrand, -np.infty, -gap, lambda x: gap,
                        lambda x: np.infty, args=(deltaE,))[0]
        part3 = dblquad(integrand, gap, np.infty, lambda x: -np.infty,
                        lambda x: -gap, args=(deltaE,))[0]
        part4 = dblquad(integrand, gap, np.infty, lambda x: gap,
                       lambda x: np.infty, args=(deltaE,))[0]

        return part1 + part2 + part3 + part4

    def integrate(deltaE):
        res = np.zeros(deltaE.shape)
        for ind,val in enumerate(deltaE):
            res[ind] = np.exp(val/T)*integrate_single(val)
        return res
    return integrate

class TunnelingRateCalculator:
    def __init__(self, deltaEmin, deltaEmax, deltaEstep, rateFunc):
        self.deltaEmin = deltaEmin
        self.deltaEmax = deltaEmax
        self.deltaEstep = deltaEstep
        self.deltaEvals = np.arange(deltaEmin, deltaEmax, deltaEstep)
        self.rateFunc = rateFunc
        self.vals = self.rateFunc(self.deltaEvals)
        self.approx = None
        self.set_approx()

    def get_approximation(self, inputs, values):
        return interp1d(inputs, values, assume_sorted=True)

    def set_approx(self):
        self.qp_approx = self.get_approximation(self.deltaEvals, self.vals)

    def increase_high_limit(self):
        new_deltaEmax = self.deltaEmax + np.abs(self.deltaEmax)
        new_inputs = np.arange(self.deltaEmax, new_deltaEmax, self.deltaEstep)
        new_vals = self.rateFunc(new_inputs)
        self.deltaEmax = new_deltaEmax
        self.deltaEvals = np.hstack((self.deltaEvals, new_inputs))
        self.vals = np.hstack((self.qp_vals, new_vals))
        self.set_approx()
        print("High limit increased, Emax= " + str(self.deltaEmax))

    def decrease_low_limit(self):
        new_deltaEmin = self.deltaEmin - np.abs(self.deltaEmin)
        new_inputs = np.arange(new_deltaEmin, self.deltaEmin, self.deltaEstep)
        new_vals = self.rateFunc(new_inputs)
        self.deltaEmin = new_deltaEmin
        self.deltaEvals = np.hstack((new_inputs, self.deltaEvals))
        self.vals = np.hstack((new_vals, self.vals))
        self.set_approx()
        print("Low limit dencreased, Emin= " + str(self.deltaEmin))

    def get_tunnling_rates(self, deltaE):
        if self.deltaEmin >= deltaE:
            self.decrease_low_limit()
        if self.deltaEmax <= deltaE:
            self.increase_high_limit()
        return self.approx(deltaE)

    def plot_rate(self):
        fig = plt.figure()
        plt.plot(self.deltaEvals, self.vals, '.')
        return fig


class NoElectronsOnDot(RuntimeError):
    pass

class DotArray:
    """
    An array of quantum dots connected to external voltage from left to
    right and to gate voltage
    """

    def __init__(self, rows, columns, VL, VR, VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, temperature, fastRelaxation=False):
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
        self.temperature = temperature
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

    def __copy__(self):
        copy_array = object.__new__(type(self))
        copy_array.fast_relaxation = self.fast_relaxation
        copy_array.rows = self.rows
        copy_array.columns = self.columns
        copy_array.VL = self.VL
        copy_array.VR = self.VR
        copy_array.VG = self.VG
        copy_array.Q = np.copy(self.Q)
        copy_array.n = np.copy(self.n)
        copy_array.CG = self.CG
        copy_array.RG = self.RG
        copy_array.Ch = self.Ch
        copy_array.Cv = self.Cv
        copy_array.Rh = self.Rh
        copy_array.Rv = self.Rv
        copy_array.R = self.R
        copy_array.temperature = self.temperature
        copy_array.totalChargePassedRight = self.totalChargePassedLeft
        copy_array.totalChargePassedLeft = self.totalChargePassedRight
        copy_array.no_tunneling_next_time = self.no_tunneling_next_time
        copy_array.invC = self.invC
        copy_array._JeigenVectors = self._JeigenVectors
        copy_array._JeigenValues = self._JeigenValues
        copy_array._JeigenVectorsInv = self._JeigenVectorsInv
        copy_array.timeStep = self.timeStep
        copy_array.default_dt = self.default_dt
        if self.fast_relaxation:
            copy_array._constQnPart = self._constQnPart
            copy_array._matrixQnPart = self._matrixQnPart
        copy_array.commonHorz = self.commonHorz
        copy_array.commonVert = self.commonVert
        copy_array.additionalLeft = self.additionalLeft
        copy_array.constWork = self.constWork
        copy_array.variableWork = np.copy(self.variableWork)
        copy_array.horizontalMatrix = self.horizontalMatrix
        copy_array.verticalMatrix = self.verticalMatrix
        copy_array._left_part_n_prime = self._left_part_n_prime
        copy_array._right_part_n_prime = self._right_part_n_prime
        copy_array.saveCurrentMap = self.saveCurrentMap
        copy_array.Ih = self.Ih
        copy_array.Iv = self.Iv
        return copy_array

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

    def setOccupation(self, n):
        self.n = n

    def setGroundCharge(self, Q):
        self.Q = Q

    def getCharge(self):
        return self.totalChargePassedRight, self.totalChargePassedLeft

    def getTimeStep(self):
        return self.timeStep

    def getCurrentMap(self):
        if self.saveCurrentMap:
            map = np.zeros((self.rows*2 - 1, self.columns+1))
            map[::2,:] = self.Ih
            if self.Iv.size:
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
        # print(-1/self._JeigenValues)
        self.timeStep = -1/np.max(self._JeigenValues)
        self.default_dt = -1/np.average(self._JeigenValues)
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
        self.commonHorz = np.pad(invCDiagMat,((0,0),(1,0)),mode='constant') + np.pad(invCDiagMat,((0,0),(0,1)),mode='constant')\
                        - lowerCDiag - upperCDiag

        lowerNCDiag = np.diag(np.copy(self.invC),k=-self.columns).reshape((self.rows-1,self.columns))
        upperNCDiag = np.diag(np.copy(self.invC),k=self.columns).reshape((self.rows-1,self.columns))
        self.commonVert = invCDiagMat[1:,:] + invCDiagMat[:-1,:] - lowerNCDiag - upperNCDiag

        self.additionalLeft = np.zeros((self.rows,self.columns+1))
        self.additionalLeft[:,0] = self.VL
        self.additionalLeft[:,-1] = -self.VR
        leftConstWork = (0.5*self.commonHorz + self.additionalLeft).flatten()

        additionalRight = -self.additionalLeft
        rightConstWork = (0.5*self.commonHorz + additionalRight).flatten()

        vertConstWork = (0.5*self.commonVert).flatten()

        self.constWork = np.hstack((rightConstWork, leftConstWork, vertConstWork, vertConstWork))
        self.variableWork = np.zeros(self.constWork.shape)
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
        variableRightWorkLen = variableRightWork.size
        variableDownWorkLen = variableDownWork.size
        self.variableWork[:variableRightWorkLen] = variableRightWork
        self.variableWork[variableRightWorkLen:2*variableRightWorkLen] = -variableRightWork
        self.variableWork[2 * variableRightWorkLen:2 * variableRightWorkLen + variableDownWorkLen] = variableDownWork
        self.variableWork[2 * variableRightWorkLen + variableDownWorkLen:] = -variableDownWork
        return self.variableWork + self.constWork

    def getRates(self):
        """
        Returns the tunnelling rate between neighboring dots
        :param T: temperature (in energy units = e^2[V])
        """
        work = self.getWork()
        if self.temperature == 0:
            work[work > 0] = 0
        else:
            notToSmall = np.abs(work) > EPS
            work[np.logical_not(notToSmall)] = -self.temperature
            work[notToSmall] = work[notToSmall]/(1 - np.exp(work[notToSmall]/self.temperature))
        return -work / self.R


    def getVoltages(self):
        return self.invC.dot(self.getNprime() + flattenToColumn(self.Q))

    def getVoltagesFromGround(self):
        return self.VG - self.Q/self.CG

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

    def tunnel(self, fromDot, toDot, charges=1):
        if fromDot[1] == self.columns:
            self.totalChargePassedRight -= charges
        elif fromDot[1] == -1:
            self.totalChargePassedLeft -= charges
        else:
            self.n[fromDot] -= 1
        if toDot[1] == self.columns:
            self.totalChargePassedRight += charges
        elif toDot[1] == -1:
            self.totalChargePassedLeft += charges
        else:
            self.n[toDot] += charges
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

    def executeAction(self, ind, charge=1):
        horzSize = self.rows*(self.columns+1)
        vertSize = (self.rows - 1)*self.columns
        if ind < horzSize: # tunnel right:
            fromDot = (ind//(self.columns+1), (ind%(self.columns+1))-1)
            toDot = (fromDot[0], fromDot[1] + 1)
            if self.saveCurrentMap:
                self.Ih[toDot] += charge
        elif ind < horzSize*2: # tunnel left
            ind -= horzSize
            fromDot = (ind//(self.columns+1), ind%(self.columns+1))
            toDot = (fromDot[0], fromDot[1] - 1)
            if self.saveCurrentMap:
                self.Ih[fromDot] -= charge
        elif ind < horzSize*2 + vertSize: # tunnel down
            ind -= horzSize*2
            fromDot = (ind//self.columns, ind%self.columns)
            toDot = (fromDot[0] + 1, fromDot[1])
            if self.saveCurrentMap:
                self.Iv[fromDot] += charge
        else: # tunnel up
            ind -= (horzSize*2 + vertSize)
            fromDot = ((ind//self.columns) + 1, ind%self.columns)
            toDot = (fromDot[0] - 1, fromDot[1])
            if self.saveCurrentMap:
                self.Iv[toDot] -= charge
        self.tunnel(fromDot, toDot)
        return fromDot, toDot

class JJArray(DotArray):
    def __init__(self, rows, columns, VL, VR, VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv,
                 temperature, scGap, fastRelaxation=False):
        DotArray.__init__(self, rows, columns, VL, VR, VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv,
                          temperature, fastRelaxation=fastRelaxation)
        self.gap = scGap
        self.Ej = self.getEj()
        self.Ec = 1/np.mean(CG)
        self.qp_rate_calculator = TunnelingRateCalculator(-1, 1, 0.01,
                                                          qp_tunneling(self.Ec, self.gap, self.temperature))
        self.cp_rate_calculator = TunnelingRateCalculator(-1, 1, 0.01,
                                                          cp_tunneling(self.Ej, self.Ec, self.temperature))
        print("Rates were calculated")
        # dbg
        # fig1 = self.qp_rate_calculator.plot_rate()
        # fig2 = self.cp_rate_calculator.plot_rate()
        # plt.show()
        # plt.close(fig1)
        # plt.close(fig2)
        # exit(0)

    def __copy__(self):
        copy_array = super().__copy__()
        copy_array.gap = self.gap
        copy_array.Ej = self.Ej
        copy_array.Ec = self.Ec
        copy_array.qp_rate_calculator = self.qp_rate_calculator
        copy_array.cp_rate_calculator = self.cp_rate_calculator
        return copy_array

    def getEj(self):
        return (self.gap/8)*np.tanh(self.gap/(2*T))

    def setConstWork(self):
        super().setConstWork()
        leftConstWork = 2*(self.commonHorz + self.additionalLeft).flatten()

        additionalRight = -self.additionalLeft
        rightConstWork = 2*(self.commonHorz + additionalRight).flatten()

        vertConstWork = (2*self.commonVert).flatten()

        self.constWork_cp = np.hstack((rightConstWork, leftConstWork, vertConstWork, vertConstWork))
        self.variableWork_cp = np.zeros(self.constWork.shape)
        self.rates = np.zeros((2*self.constWork.shape[0],))

    def getWork(self):
        qp_work = super().getWork()
        cp_work = 2*self.variableWork + self.variableWork_cp
        return qp_work, cp_work

    def getRates(self):
        """
        Returns the tunnelling rate between neighboring dots
        :param T: temperature (in energy units = e^2[V])
        """
        qp_work, cp_work = self.getWork()
        if self.temperature == 0:
            raise NotImplementedError
        else:
            mid = self.rates.size // 2
            self.rates[:mid] = self.qp_rate_calculator.get_tunnling_rates(-qp_work) / self.R
            self.rates[mid:] = self.cp_rate_calculator.get_tunnling_rates(-cp_work) / (self.R**2)
        return self.rates

    def executeAction(self, ind, charge=1):
        horzSize = self.rows*(self.columns+1)
        vertSize = (self.rows - 1)*self.columns
        if ind < 2*(horzSize + vertSize):
            fromDot, toDot = super().executeAction(ind, charge=1)
        else:
            ind = ind - 2*(horzSize + vertSize)
            fromDot, toDot = super().executeAction(ind, charge=2)
        return fromDot, toDot

class Simulator:
    """
    Preforms Gillespie simulation on an array of quantum dots
    :param rows: number of rows in the array
    :param columns: number of columns in the array
    """
    def __init__(self, index, VL0, VR0, Q0, n0, dotArray):
        self.dotArray = copy(dotArray)
        self.VL = VL0
        self.VR = VR0
        self.index = index
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

    def calcCurrent(self, t,print_stats=False,fullOutput=False, currentMap=False):
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
            dt = self.executeStep(printState=print_stats)
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
        ns = []
        Qs = []
        t = self.dotArray.getTimeStep()/10
        for r in range(rep):
            # running once to get to steady state
            # self.calcCurrent(t)
            # now we are in steady state calculate current
            rightCurrent, leftCurrnet, n, Q = self.calcCurrent(t, fullOutput=True)
            Ileft.append(leftCurrnet)
            Iright.append(rightCurrent)
            ns.append(n)
            Qs.append(Q)
        ns = np.array(ns)
        Qs = np.array(Qs)
        x = t*np.arange(rep)
        print("Mean right: " + str(np.average(Iright))+ " Var right: " + str(np.var(Iright)))
        print("Mean left: " + str(np.average(Ileft))+ " Var left: " + str(np.var(Ileft)))
        print("Tail mean right: " + str(np.average(Iright[len(Iright)//4:])) + " Tail var right: " + str(np.var(Iright[len(Iright)//4:])))
        print("Tail mean left: " + str(np.average(Ileft[len(Ileft)//4:])) + " Tail var left: " + str(np.var(Ileft[len(Ileft)//4:])))
        plt.figure()
        plt.plot(x,Ileft, x,Iright, x, (np.array(Ileft) - np.array(Iright))**2)
        plt.figure()
        plt.plot((np.array(Ileft) - np.array(Iright)) ** 2)
        plt.figure()
        plt.plot(x, ns[:,0,:], x, ns[:,1,:])
        plt.figure()
        plt.plot(x, Qs[:,0,:], x, Qs[:,1,:])
        plt.show()

    def saveState(self, I, V, n=None, Q=None, Imaps=None, fullOutput=False, currentMap=False, basePath=''):
        baseName = basePath + "_temp_" + str(self.index)
        np.save(baseName + "_I", np.array(I))
        np.save(baseName + "_Vind", np.array(V))
        if fullOutput:
            np.save(baseName + "_ns", np.array(n))
            np.save(baseName + "_Qs", np.array(Q))
        np.save(baseName + "_n", self.dotArray.getOccupation())
        np.save(baseName + "_Q", self.dotArray.getGroundCharge())
        if currentMap:
            np.save(baseName + "_current_map", np.array(Imaps))

    def loadState(self,  fullOutput=False, currentMap=False, basePath=''):
        baseName = basePath + "_temp_" + str(self.index)
        I = np.load(baseName + "_I.npy")
        Vind = np.load(baseName + "_Vind.npy")
        n = np.load(baseName + "_n.npy")
        Q = np.load(baseName + "_Q.npy")
        res = (I,Vind,n,Q)
        if fullOutput:
            ns = np.load(baseName + "_ns.npy")
            Qs = np.load(baseName + "_Qs.npy")
            res = res + (ns, Qs)
        if currentMap:
            Imaps = np.load(baseName + "_current_map.npy")
            res = res + (Imaps,)
        return res

    def clearState(self, fullOutput=False, currentMap=False, basePath=''):
        baseName = basePath + "_temp_" + str(self.index)
        os.remove(baseName + "_I.npy")
        os.remove(baseName + "_Vind.npy")
        os.remove(baseName + "_n.npy")
        os.remove(baseName + "_Q.npy")
        if fullOutput:
            os.remove(baseName + "_ns.npy")
            os.remove(baseName + "_Qs.npy")
        if currentMap:
            os.remove(baseName + "_current_map.npy")
        return True

    def calcIV(self, Vmax, Vstep, vSym, fullOutput=False, print_stats=False, currentMap=False, basePath="", resume=False):
        I = []
        ns = []
        Qs = []
        Imaps = []
        tStep = self.dotArray.getTimeStep()
        if vSym:
            Vstep /= 2
            Vmax /= 2
            VR_vec = np.arange(self.VR-(self.VL/2), self.VR - Vmax, -Vstep)
            VR_vec = np.hstack((VR_vec, np.flip(VR_vec)))
            VL_vec = np.arange(self.VL/2+self.VR, Vmax + self.VR, Vstep)
            VL_vec = np.hstack((VL_vec, np.flip(VL_vec)))
        else:
            VL_vec = np.arange(self.VL, Vmax+self.VR, Vstep)
            VL_vec = np.hstack((VL_vec, np.flip(VL_vec)))
            VR_vec = self.VR * np.ones(VL_vec.shape)
        VL_res = np.copy(VL_vec)
        VR_res = np.copy(VR_vec)
        Vind_addition = 0
        if resume:
            resumeParams = self.loadState(fullOutput=fullOutput, currentMap=currentMap, basePath=basePath)
            I = list(resumeParams[0])
            Vind = int(resumeParams[1])
            VL_vec = VL_vec[Vind+1:]
            VR_vec = VR_vec[Vind+1:]
            Vind_addition = Vind + 1
            self.dotArray.setOccupation(resumeParams[2])
            self.dotArray.setGroundCharge(resumeParams[3])
            if fullOutput:
                ns = list(resumeParams[4])
                Qs = list(resumeParams[5])
            if currentMap:
                Imaps = list(resumeParams[6])
        Vind = 0
        for VL,VR in zip(VL_vec, VR_vec):
            self.dotArray.changeVext(VL, VR)
            # running once to get to steady state
            self.calcCurrent(tStep/10)
            # now we are in steady state calculate current
            stepRes = self.calcCurrent(tStep, print_stats=print_stats, fullOutput=fullOutput, currentMap=currentMap)
            rightCurrent = stepRes[0]
            leftCurrnet = stepRes[1]
            if fullOutput:
                ns.append(stepRes[2])
                Qs.append(stepRes[3])
            if currentMap:
                Imaps.append(stepRes[-1])
            I.append((rightCurrent+leftCurrnet)/2)
            print(VL - VR, end=',')
            self.saveState(I, [Vind + Vind_addition], ns, Qs, Imaps, fullOutput=fullOutput,
                           currentMap=currentMap, basePath=basePath)
        res = (np.array(I), VL_res - VR_res)
        if fullOutput:
            res = res + (np.array(ns), np.array(Qs))
        if currentMap:
            res = res + (np.array(Imaps),)
        self.clearState(fullOutput=fullOutput, currentMap=currentMap, basePath=basePath)
        return res

    def printState(self):
        self.dotArray.printState()

class GraphSimulator:
    """
    Calculating steady state current for an array of quantum dots by using linear programming
    """
    def __init__(self, index, VL0, VR0, Q0, n0, dotArray):
        self.dotArray = copy(dotArray)
        self.n0 = n0
        self.QG = Q0
        self.VL = VL0
        self.VR = VR0
        self.counter = 0
        self.index = index

        self.edgesMat = None
        self.states = None
        self.prob = None
        self.lyaponuv = None
        self.rates_diff_left = None
        self.rates_diff_right = None

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
        edges_line = [0]
        while current_state_ind < len(states):
            self.dotArray.n = self.reshape_to_array(np.copy(states[current_state_ind]))
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
            edges_line = [0] * len(edges_line)
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
        second_round = False
        sucess = False
        while not sucess:
            for i in range(null.shape[1]):
                candidate = null[:,i]
                if second_round:
                    candidate[np.abs(candidate) < EPS] = 0
                if (candidate >= 0).all() or (candidate <= 0).all():
                    sol.append(candidate / np.sum(candidate))
            if len(sol) > 0:
                self.prob = sol[0]
                sucess = True
            elif second_round:
                self.prob = np.zeros((len(self.states),))
                self.prob[0] = 1
                sucess = True
            else:
                second_round = True
        return True

    def get_average_state(self, Q):
        self.find_probabilities(Q)
        average_state = np.sum(np.multiply(self.states.T, self.prob),axis=1)
        return self.reshape_to_array(average_state)

    def get_average_voltages(self, Q):
        n_copy = np.copy(self.dotArray.n)
        Q_copy = np.copy(self.dotArray.Q)
        n = self.get_average_state(Q)
        self.dotArray.n = self.reshape_to_array(n)
        self.dotArray.Q = self.reshape_to_array(Q)
        v = self.dotArray.getVoltages()
        self.dotArray.n = n_copy
        self.dotArray.Q = Q_copy
        return v

    def get_voltages_from_ground(self, Q):
        Q_copy = np.copy(self.dotArray.Q)
        self.dotArray.Q = self.reshape_to_array(Q)
        v = self.dotArray.getVoltagesFromGround()
        self.dotArray.Q = Q_copy
        return v

    def calc_lyaponuv(self, Q):
        n_copy = np.copy(self.dotArray.n)
        Q_copy = np.copy(self.dotArray.Q)
        n = self.get_average_state(Q)
        self.dotArray.n = self.reshape_to_array(n)
        self.dotArray.Q = self.reshape_to_array(Q)
        diff_from_equi = np.sum(flattenToColumn(Q) - self.dotArray.get_steady_Q_for_n())
        self.dotArray.n = n_copy
        self.dotArray.Q = Q_copy
        return diff_from_equi

    def calc_on_grid(self, Qmin, Qmax, f, res_len=1):
        Qmin = Qmin.flatten()
        Qmax = Qmax.flatten()
        dq = DQ
        coordinates = [np.arange(Qmin[i], Qmax[i], dq) for i in range(Qmin.size)]
        grid = np.meshgrid(*coordinates, indexing='ij')
        grid_array = np.moveaxis(np.array(grid),0,-1)
        res = [np.zeros(grid[0].shape) for i in range(res_len)]
        it = np.nditer(grid[0], flags=['multi_index'])
        while not it.finished:
            index = it.multi_index
            curr_Q = grid_array[index]
            for i in range(res_len):
                f_res = f(curr_Q).flatten()
                res[i][index] = f_res[i]
            it.iternext()
        return grid_array, res

    def set_lyaponuv(self, Qmin, Qmax):
        self.Q_grid, diff_from_eq = self.calc_on_grid(Qmin, Qmax, self.calc_lyaponuv)
        res = diff_from_eq[0]
        for axis in range(len(res.shape)):
            res = cumtrapz(res, dx=DQ, axis=axis, initial=0)
        self.lyaponuv = res

    def calc_lyaponuv_grad(self, Q):
        n = self.get_average_state(self.reshape_to_array(Q))
        self.dotArray.n = n
        self.dotArray.Q = self.reshape_to_array(Q)
        return Q.flatten() - self.dotArray.get_steady_Q_for_n().flatten()

    def plot_average_voltages(self, Qmin, Qmax):
        from mayavi import mlab
        Q_grid, voltages = self.calc_on_grid(Qmin, Qmax, self.get_average_voltages, res_len=Qmin.size)
        _, voltages_from_ground = self.calc_on_grid(Qmin, Qmax, self.get_voltages_from_ground, res_len=Qmin.size)
        fig1 = mlab.figure()
        surf1 = mlab.surf(Q_grid[:, :, 0], Q_grid[:, :, 1],voltages[0], colormap='Blues')
        surf2 = mlab.surf(Q_grid[:, :, 0], Q_grid[:, :, 1],voltages_from_ground[0], colormap='Oranges')
        fig1 = mlab.figure()
        surf1 = mlab.surf(Q_grid[:, :, 0], Q_grid[:, :, 1], voltages[1], colormap='Blues')
        surf2 = mlab.surf(Q_grid[:, :, 0], Q_grid[:, :, 1], voltages_from_ground[1], colormap='Oranges')
        mlab.show()

    def find_next_QG_using_lyaponuv(self, basePath):
        # q_shift = Q_SHIFT
        # peaks = (np.array([]),)
        # while not peaks[0].size:
        #     self.set_lyaponuv(self.QG - q_shift, self.QG + q_shift)
        #     q_shift *= 2
        #     peaks = detect_local_minima(self.lyaponuv)
        # Qind = np.argmin(np.sum((self.Q_grid[peaks] - self.QG)**2,axis=1))
        # self.QG = self.Q_grid[peaks][Qind]
        # dbg
        self.set_lyaponuv(self.QG - Q_SHIFT, self.QG + Q_SHIFT)
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(self.Q_grid[:,:,0],self.Q_grid[:,:,1],self.lyaponuv)
        x = np.searchsorted(self.Q_grid[:,0,0], self.QG[0,0])
        y = np.searchsorted(self.Q_grid[0,:,1], self.QG[0,1])
        ax.scatter3D(self.QG[0,0], self.QG[0,1], self.lyaponuv[x,y], marker='o',color='red')
        str_ind = "0"*(4-len(str(self.counter))) + str(self.counter)
        plt.savefig(basePath + "Lyaponuv_" + str_ind)
        plt.close(fig)
        self.counter += 1

    def find_next_QG_using_gradient_descent(self):
        flag = False
        rep = 0
        lr = INI_LR
        while (not flag) and rep < GRAD_REP:
            res, flag = simple_gadient_descent(self.calc_lyaponuv_grad, self.QG, lr=lr, plot_lc=False)
            rep += 1
            lr = lr/10
        if flag: #if gradient descent did not converge skipping the point
            self.QG = self.reshape_to_array(res)

    def calcCurrent(self, fullOutput=False, basePath=""):
        self.find_next_QG_using_gradient_descent()
        #dbg
        self.find_next_QG_using_lyaponuv(basePath)
        print(self.QG)
        self.plot_average_voltages(self.QG-Q_SHIFT, self.QG+Q_SHIFT)
        #dbg
        n_avg = self.reshape_to_array(self.get_average_state(self.QG))
        self.n0 = np.floor(n_avg)
        left_current = np.sum(self.prob*self.rates_diff_left.T)
        right_current = np.sum(self.prob*self.rates_diff_right.T)
        if fullOutput:
            return right_current, left_current, n_avg, self.reshape_to_array(self.QG)
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

    def saveState(self, I, Vind, n=None, Q=None, fullOutput=False, basePath=''):
        baseName = basePath + "_temp_" + str(self.index)
        np.save(baseName + "_I", np.array(I))
        np.save(baseName + "_Vind", np.array(Vind))
        if fullOutput:
            np.save(baseName + "_ns", np.array(n))
            np.save(baseName + "_Qs", np.array(Q))
        np.save(baseName + "_n", self.n0)
        np.save(baseName + "_Q", self.QG)

    def loadState(self, fullOutput=False, basePath=''):
        baseName = basePath + "_temp_" + str(self.index)
        I = np.load(baseName + "_I.npy")
        Vind = np.load(baseName + "_Vind.npy")
        n = np.load(baseName + "_n.npy")
        Q = np.load(baseName + "_Q.npy")
        res = (I,Vind, n, Q)
        if fullOutput:
            ns = np.load(baseName + "_ns.npy")
            Qs = np.load(baseName + "_Qs.npy")
            res = res + (ns, Qs)
        return res

    def removeState(self, fullOutput=False, basePath=''):
        baseName = basePath + "_temp_" + str(self.index)
        os.remove(baseName + "_I.npy")
        os.remove(baseName + "_Vind.npy")
        os.remove(baseName + "_n.npy")
        os.remove(baseName + "_Q.npy")
        if fullOutput:
            os.remove(baseName + "_ns.npy")
            os.remove(baseName + "_Qs.npy")
        return True

    def calcIV(self, Vmax, Vstep, vSym, fullOutput=False, print_stats=False, currentMap=False,
               basePath="", resume=False):
        I = []
        ns = []
        Qs = []
        Vind_addition = 0
        if vSym:
            Vstep /= 2
            Vmax /= 2
            VR_vec = np.arange(self.VR, -Vmax + self.VR, -Vstep)
            VR_vec = np.hstack((VR_vec, np.flip(VR_vec)))
        VL_vec = np.arange(self.VL, Vmax + self.VR, Vstep)
        VL_vec = np.hstack((VL_vec, np.flip(VL_vec)))
        VL_res = np.copy(VL_vec)
        if not vSym:
            VR_vec = self.VR*np.ones(VL_vec.shape)
        VR_res = np.copy(VR_vec)
        if resume:
            resumeParams = self.loadState(fullOutput=fullOutput, basePath=basePath)
            I = list(resumeParams[0])
            Vind = int(resumeParams[1])
            VL_vec = VL_vec[Vind+1:]
            VR_vec = VR_vec[Vind+1:]
            Vind_addition = Vind + 1
            self.n0 = resumeParams[2]
            self.QG = resumeParams[3]
            if fullOutput:
                ns = list(resumeParams[4])
                Qs = list(resumeParams[5])
        Vind = 0
        for VL,VR in zip(VL_vec,VR_res):
            print(VL-VR,end=',')
            self.dotArray.changeVext(VL, VR)
            res = self.calcCurrent(fullOutput=fullOutput, basePath=basePath)
            rightCurrent = res[0]
            leftCurrent = res[1]
            if fullOutput:
                n = res[2]
                Q = res[3]
                ns.append(n)
                Qs.append(Q)
            I.append((rightCurrent + leftCurrent) / 2)
            self.saveState(I, [Vind + Vind_addition], ns, Qs, fullOutput=fullOutput, basePath=basePath)
            Vind += 1
        result = (np.array(I), VL_res - VR_res)
        if fullOutput:
            result = result + (ns, Qs)
        self.removeState(fullOutput=fullOutput, basePath=basePath)
        return result

def runSingleSimulation(index, VL0, VR0, vSym, Q0, n0,Vmax, Vstep, dotArray,
                        fullOutput=False, printState=False, useGraph=False, currentMap=False,
                        basePath="", resume=False):
    if useGraph:
        simulator = GraphSimulator(index, VL0, VR0, Q0, n0, dotArray)
    else:
        simulator = Simulator(index, VL0, VR0, Q0, n0, dotArray)
    out = simulator.calcIV(Vmax, Vstep, vSym, fullOutput=fullOutput, print_stats=printState,
                           currentMap=currentMap, basePath=basePath, resume=resume)
    array_params = simulator.getArrayParameters()
    return out + (array_params,)

def saveRandomParams(VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv,basePath):
    baseName = basePath + "_temp_"
    np.save(baseName + "VG", VG)
    np.save(baseName + "Q0", Q0)
    np.save(baseName + "n0", n0)
    np.save(baseName + "CG", CG)
    np.save(baseName + "RG", RG)
    np.save(baseName + "Ch", Ch)
    np.save(baseName + "Cv", Cv)
    np.save(baseName + "Rh", Rh)
    np.save(baseName + "Rv", Rv)
    return True

def loadRandomParams(basePath):
    baseName = basePath + "_temp_"
    VG0 = np.load(baseName + "VG.npy")
    Q0 = np.load(baseName + "Q0.npy")
    n0 = np.load(baseName + "n0.npy")
    CG = np.load(baseName + "CG.npy")
    RG = np.load(baseName + "RG.npy")
    Ch = np.load(baseName + "Ch.npy")
    Cv = np.load(baseName + "Cv.npy")
    Rh = np.load(baseName + "Rh.npy")
    Rv = np.load(baseName + "Rv.npy")
    return VG0, Q0, n0, CG, RG, Ch, Cv, Rh, Rv

def removeRandomParams(basePath):
    baseName = basePath + "_temp_"
    os.remove(baseName + "VG.npy")
    os.remove(baseName + "Q0.npy")
    os.remove(baseName + "n0.npy")
    os.remove(baseName + "CG.npy")
    os.remove(baseName + "RG.npy")
    os.remove(baseName + "Ch.npy")
    os.remove(baseName + "Cv.npy")
    os.remove(baseName + "Rh.npy")
    os.remove(baseName + "Rv.npy")
    return True

def runFullSimulation(VL0, VR0, vSym, VG0, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, rows, columns,
                      Vmax, Vstep, temperature=0, repeats=1, savePath=".", fileName="", fullOutput=False,
                      printState=False, checkSteadyState=False, useGraph=False, fastRelaxation=False,
                      currentMap=False, dbg=False, plotCurrentMaps=False, resume=False, superconducting=False,
                      gap=0):
    basePath = os.path.join(savePath, fileName)
    if plotCurrentMaps:
        print("Plotting Current Maps")
        avgImaps = np.load(basePath + "_Imap.npy")
        V = np.load(basePath + "_V.npy")
        saveCurrentMaps(avgImaps, V, basePath + "_Imap")
        exit(0)
    if not resume:
        print("Saving array parameters")
        saveRandomParams(np.array(VG0),
                         np.array(Q0), np.array(n0), np.array(CG),
                         np.array(RG), np.array(Ch), np.array(Cv),
                         np.array(Rh), np.array(Rv), basePath)
    else:
        print("Loading array parameters")
        VG0, Q0, n0, CG, RG, Ch, Cv, Rh, Rv = loadRandomParams(basePath)
    if superconducting:
        prototypeArray = JJArray(rows, columns, VL0, VR0, np.array(VG0), np.array(Q0), np.array(n0), np.array(CG),
                                 np.array(RG), np.array(Ch), np.array(Cv), np.array(Rh), np.array(Rv),
                                 temperature, gap, fastRelaxation=fastRelaxation)
        print("Superconducting prototype array was created")
    else:
        prototypeArray = DotArray(rows, columns, VL0, VR0, np.array(VG0), np.array(Q0), np.array(n0), np.array(CG),
                                  np.array(RG), np.array(Ch), np.array(Cv), np.array(Rh), np.array(Rv),
                                  temperature, fastRelaxation=fastRelaxation)
        print("Normal prototype array was created")
    if checkSteadyState:
        print("Running steady state convergance check")
        simulator = Simulator(0, VL0, VR0, Q0, n0, prototypeArray)
        simulator.checkSteadyState(repeats)
        simulator.dotArray.changeVext(VL0+Vstep, VR0)
        simulator.checkSteadyState(repeats)
        exit(0)
    Is = []
    ns = []
    Qs = []
    Imaps = []
    if not dbg:
        print("Starting parallel run")
        pool = Pool(processes=repeats)
        results = []
        for repeat in range(repeats):
            res = pool.apply_async(runSingleSimulation,
                                    (repeat, VL0, VR0, vSym, Q0, n0, Vmax, Vstep, prototypeArray, fullOutput,
                                     printState, useGraph, currentMap,basePath, resume))
            results.append(res)
        for res in results:
            result = res.get()
            I = result[0]
            print(I)
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

    else: #dbg
        print("Starting serial run")
        for repeat in range(repeats):
            result = runSingleSimulation(repeat, VL0, VR0, vSym, Q0, n0, Vmax, Vstep, prototypeArray, fullOutput,
                                         printState, useGraph, currentMap,basePath, resume)
            I = result[0]
            currentMapInd = 2
            if fullOutput:
                n = result[2]
                Q = result[3]
                ns.append(n)
                Qs.append(Q)
                currentMapInd = 4
            if currentMap:
                Imaps.append(result[currentMapInd])
            Is.append(I)
        V = result[1]
        params = result[-1]
        # dbg
    print("Saving results")
    avgI = np.mean(np.array(Is), axis=0)
    if fullOutput:
        avgN = np.mean(np.array(ns), axis=0)
        avgQ = np.mean(np.array(Qs), axis=0)
        np.save(basePath + "_n", avgN)
        np.save(basePath + "_Q", avgQ)
        np.save(basePath + "_full_I", np.array(Is))
    fig = plt.figure()
    plt.plot(V[:V.size//2], avgI[:V.size//2], '.b',V[V.size//2:],
             avgI[V.size//2:],
             '.r')
    plt.savefig(basePath + "_IV.png")
    np.save(basePath + "_I", avgI)
    np.save(basePath + "_V", V)
    plt.close(fig)
    if currentMap:
        avgImaps = np.mean(np.array(Imaps),axis=0)
        np.save(basePath + "_Imap", avgImaps)
    removeRandomParams(basePath)
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
    parser.add_option("-T", "--temperature", dest="T",
                      help="Environment temperature (in units of planckConstant/timeUnits) [default: %default]",
                      default=0, type=float)
    parser.add_option("--gap", dest="gap",
                      help="superconducting gap (in units of planckConstant/timeUnits)[default: %default]",
                      default=0, type=float)
    parser.add_option("-M", "--height", dest="M", help="number of lines in "
                      "the array [default: %default]", default=1, type=int)
    parser.add_option("-N", "--width", dest="N", help="number of columns in "
                      "the array [default: %default]", default=1, type=int)
    parser.add_option("--vr", dest="VR", help="right electrode voltage (in units of"
                                              " planckConstant/electronCharge*timeUnits) [default: %default]",
                      default=0, type=float)
    parser.add_option("--vmin", dest="Vmin", help="minimum external voltage  (in units of"
                                              " planckConstant/electronCharge*timeUnits)"
                      " [default: %default]", default=0, type=float)
    parser.add_option("--vmax", dest="Vmax", help="maximum external voltage  (in units of"
                                              " planckConstant/electronCharge*timeUnits)"
                      " [default: %default]", default=10, type=float)
    parser.add_option("--vstep", dest="vStep", help="size of voltage step  (in units of"
                                              " planckConstant/electronCharge*timeUnits)[default: %default]",
                      default=1, type=float)
    parser.add_option("--symmetric-v", dest="vSym", help="Voltage raises symmetric on VR and VL["
                                                    "default: %default]", default=False, action='store_true')
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
    parser.add_option("--current-map", dest="current_map", help="if true fraes for clip of current distribution during"
                                                                " simulation will be created and saved [Default:%default]",
                      default=False, action='store_true')
    parser.add_option("--plot-current-map", dest="plot_current_map", help="if true clip of current distribution will"
                                                                          " be plotted using former saved frames (from"
                                                                          " a former run with same file name and"
                                                                          " location and the flag --current-map"
                                                                          " [Default:%default]",
                      default=False, action='store_true')
    parser.add_option("--dbg", dest="dbg", help="Avoids parallel running for debugging [Default:%default]",
                      default=False, action='store_true')
    parser.add_option("--resume", dest="resume", help="Resume failed run from last checkpoint [Default:%default]",
                      default=False, action='store_true')
    parser.add_option("--superconducting", dest="sc", help="use superconducting array [Default:%default]",
                      default=False, action='store_true')
    parser.add_option("-o", "--output-folder", dest="output_folder",
                      help="Output folder [default: current folder]",
                      default='.')

    # Disorder Parameters
    parser.add_option("--vg-avg", dest="VG_avg", help="Gate voltage average  (in units of"
                                                      " planckConstant/electronCharge*timeUnits) [default: %default]",
                      default=1, type=float)
    parser.add_option("--vg-std", dest="VG_std", help="Gate voltage std  (in units of"
                                              " planckConstant/electronCharge*timeUnits) [default: %default]",
                      default=0, type=float)
    parser.add_option("--c-avg", dest="C_avg", help="capacitance of junctions average (in units of"
                                                    " timeUnits*electronCharge^2/planckConstant) [default: %default]",
                      default=1, type=float)
    parser.add_option("--c-std", dest="C_std", help="capacitance of junctions std (in units of"
                                                    " timeUnits*electronCharge^2/planckConstant) [default: %default]",
                      default=0, type=float)
    parser.add_option("--cg-avg", dest="CG_avg", help="Gate Capacitors capacitance average (in units of"
                                                    " timeUnits*electronCharge^2/planckConstant) [default: %default]",
                      default=1, type=float)
    parser.add_option("--cg-std", dest="CG_std", help="Gate Capacitors capacitance std (in units of"
                                                    " timeUnits*electronCharge^2/planckConstant) [default: %default]",
                      default=0, type=float)
    parser.add_option("--r-avg", dest="R_avg", help="junctions resistance average (in units of"
                                                    " planckConstant/electronCharge^2) [default: %default]",
                      default=1, type=float)
    parser.add_option("--r-std", dest="R_std", help="junctions resistance std (in units of"
                                                    " planckConstant/electronCharge^2) [default: %default]",
                      default=0, type=float)
    parser.add_option("--custom-rh", dest="custom_rh", help="list of r horizontal values ordered as numpy array."
                                                            " Overrides random r parameters [default: %default]",
                      default="")
    parser.add_option("--custom-rv", dest="custom_rv", help="list of r vertical values ordered as numpy array."
                                                             " Overrides random r parameters [default: %default]",
                      default="")
    parser.add_option("--custom-ch", dest="custom_ch", help="list of c horizontal values ordered as numpy array."
                                                            " Overrides random c parameters [default: %default]",
                      default="")
    parser.add_option("--custom-cv", dest="custom_cv", help="list of c vertical values ordered as numpy array."
                                                            " Overrides random c parameters [default: %default]",
                      default="")
    parser.add_option("--rg-avg", dest="RG_avg", help="Gate Resistors resistance average (in units of"
                                                    " planckConstant/electronCharge^2) [default: %default]",
                      default=1, type=float)
    parser.add_option("--rg-std", dest="RG_std", help="Gate Resistors resistance std (in units of"
                                                    " planckConstant/electronCharge^2) [default: %default]",
                      default=0, type=float)
    parser.add_option("--n-avg", dest="n0_avg", help="initial number of "
                      "electrons on each dot average [default:%default]",
                      default=0, type=float)
    parser.add_option("--n-std", dest="n0_std", help="initial number of "
                      "electrons on each dot std [default:%default]",
                      default=0, type=float)
    parser.add_option("--q-avg", dest="Q0_avg", help="initial charge on gate capacitors average (in units of"
                                                    " electronCharge) [default:%default]",
                      default=0, type=float)
    parser.add_option("--q-std", dest="Q0_std", help="initial charge on gate capacitors std (in units of"
                                                    " electronCharge) [default:%default]",
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
    elif dist == 'uniform' or (dist == 'exp' and not only_positive):
        res = np.random.uniform(low=avg-std, high=avg+std,size=(M,N))
    elif dist == 'two_points':
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
    elif dist == 'exp':
        r = np.random.uniform(low=np.log2(max(avg - std,0.01)), high=np.log2(max(avg + std,0.01)), size=(M, N))
        res = 2**r
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
    T = options.T
    gap = options.gap
    sc = options.sc
    VG = create_random_array(rows, columns, options.VG_avg, options.VG_std, dist,
                             False)
    Q0 = create_random_array(rows, columns, options.Q0_avg, options.Q0_std, dist,
                             False)
    n0 = create_random_array(rows, columns, options.n0_avg, options.n0_std, dist,False)
    CG = create_random_array(rows, columns, options.CG_avg, options.CG_std, dist,True)
    RG = create_random_array(rows, columns, options.RG_avg, options.RG_std, dist,True)
    if options.custom_ch.replace('\"', '') and options.custom_cv.replace('\"', ''):
        Ch = literal_eval(options.custom_ch.replace('\"', ''))
        Cv = literal_eval(options.custom_cv.replace('\"', ''))
    else:
        Ch = create_random_array(rows, columns + 1, options.C_avg, options.C_std, dist,
                                 True)
        Cv = create_random_array(rows - 1, columns, options.C_avg, options.C_std, dist,
                                 True)
    if options.custom_rh.replace('\"','') and options.custom_rv.replace('\"',''):
        Rh = literal_eval(options.custom_rh.replace('\"',''))
        Rv = literal_eval(options.custom_rv.replace('\"',''))
    else:
        Rh = create_random_array(rows, columns + 1, options.R_avg, options.R_std, dist,
                            True)
        Rv = create_random_array(rows - 1, columns, options.R_avg, options.R_std, dist,
                            True)
    Vmax = options.Vmax
    Vstep = options.vStep
    vSym = options.vSym
    repeats = options.repeats
    savePath = options.output_folder
    fileName = options.fileName
    fullOutput = options.fullOutput
    use_graph = options.use_graph
    fast_relaxation = options.RG_avg * options.CG_avg < options.C_avg * options.R_avg
    current_map = options.current_map
    dbg = options.dbg
    resume = options.resume
    plot_current_map = options.plot_current_map
    # Running Simulation
    if not os.path.exists(savePath):
        os.mkdir(savePath)
    elif not os.path.isdir(savePath):
        print("the given path exists but is a file")
        exit(0)
    import cProfile, pstats, io
    from pstats import SortKey
    pr = cProfile.Profile()
    pr.enable()
    array_params = runFullSimulation(VL0, VR0, vSym, VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, rows,  columns,
                                     Vmax, Vstep, temperature=T, repeats=repeats, savePath=savePath, fileName=fileName,
                                     fullOutput=fullOutput, printState=False, useGraph=use_graph,
                                     fastRelaxation=fast_relaxation, currentMap=current_map,
                                     dbg=dbg, plotCurrentMaps=plot_current_map, resume=resume,
                                     checkSteadyState=False, superconducting=sc, gap=gap)
    pr.disable()
    s = io.StringIO()
    sortby = SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print(s.getvalue())
    saveParameters(savePath, fileName, options, array_params)
    exit(0)
