__author__ = 'shahar'

import numpy as np
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
from matplotlib import pyplot as plt
import matplotlib.animation as animation
# from mayavi import mlab
from optparse import OptionParser
import re
import os
from multiprocessing import Pool
from scipy.linalg import null_space
from scipy.integrate import cumtrapz, dblquad
from scipy.interpolate import interp1d
from scipy.stats import norm
from mpl_toolkits.mplot3d import Axes3D
from ast import literal_eval
from copy import copy


EPS = 0.0001
# Gillespie Constants
MIN_STEPS = 1000
STEADY_STATE_VAR = 1e-4
ALLOWED_ERR = 1e-4
STEADY_STATE_REP = 10
INV_DOS = 0.01

# Tau Leaping Constants
TAU_EPS = 0.03
# Lyaponuv Constants
DQ = 0.1
Q_SHIFT = 10
GRAD_REP = 5
INI_LR = 0.1


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
def simple_gadient_descent(grad, x0, eps=1e-4, lr=1e-2, max_iter=100, plot_lc=True):
    x=x0.flatten()
    curr_grad = grad(x)
    if plot_lc:
        gradvec = []
    iter=0
    while np.max(np.abs(curr_grad)) > eps and iter < max_iter:
        if plot_lc:
            print(iter)
            gradvec.append(curr_grad)
        x = x - curr_grad*lr
        curr_grad = grad(x)
        iter += 1
    if plot_lc:
        plt.figure()
        plt.plot(gradvec)
        plt.show()
    return x, iter < max_iter

# Methods for calculating superconducting related tunneling rates

def high_impadance_p(x,Ec,T,kappa):
    sigma = 4*(kappa**2)*Ec*T
    mu = kappa**2*Ec
    return norm.pdf(x, loc=mu, scale=np.sqrt(sigma))

def fermi_dirac_dist(x,T):
    return 1/(1 + np.exp(x/T))

def qp_density_of_states(x,energy_gap):
    return np.abs(x) / np.sqrt(x**2-energy_gap**2)

def cp_tunneling(x, Ec, Ej, T):
    return (np.pi*Ej)**2 * high_impadance_p(x, Ec, T, 2)

def qp_integrand(x2, x1, deltaE, Ec, gap, T):
    return qp_density_of_states(x1, gap)*qp_density_of_states(x2, gap)*fermi_dirac_dist(-x1, T)\
           *fermi_dirac_dist(x2, T)*high_impadance_p(x2-x1 - deltaE, Ec, T, 1)

def qp_tunneling_single(deltaE, Ec, gap, T):
    part1 = dblquad(qp_integrand, -np.infty, -gap, lambda x: -np.infty,
                    lambda x: -gap, args=(deltaE,Ec, gap, T))[0]
    part2 = dblquad(qp_integrand, -np.infty, -gap, lambda x: gap,
                    lambda x: np.infty, args=(deltaE,Ec, gap, T))[0]
    part3 = dblquad(qp_integrand, gap, np.infty, lambda x: -np.infty,
                    lambda x: -gap, args=(deltaE,Ec, gap, T))[0]
    part4 = dblquad(qp_integrand, gap, np.infty, lambda x: gap,
                    lambda x: np.infty, args=(deltaE,Ec, gap, T))[0]
    return part1 + part2 + part3 + part4

def qp_tunneling(deltaE, Ec, gap, T):
        res = np.zeros(deltaE.shape)
        for ind,val in enumerate(deltaE):
            res[ind] = np.exp(val/T)*qp_tunneling_single(val, Ec, gap, T)
        return res

class TunnelingRateCalculator:
    def __init__(self, deltaEmin, deltaEmax, deltaEstep, rateFunc, Ec, T, otherParam, dirPath):
        self.T = T
        self.otherParam = otherParam
        self.Ec = Ec
        self.dirName = os.path.join(dirPath, "tunneling_rate_" + str(self.T) + "_" + str(self.otherParam))
        if not os.path.isdir(dirPath):
            os.mkdir(dirPath)
        self.deltaEmin = deltaEmin
        self.deltaEmax = deltaEmax
        self.deltaEstep = deltaEstep
        self.rateFunc = rateFunc
        self.set_results()
        self.set_approx()

    def set_results(self):
        can_load = False
        deltaEmin = None
        deltaEmax = None
        deltaEstep = None
        if os.path.isdir(self.dirName):
            deltaEmin = np.load(os.path.join(self.dirName, "deltaEmin.npy"))
            deltaEmax = np.load(os.path.join(self.dirName, "deltaEmax.npy"))
            deltaEstep = np.load(os.path.join(self.dirName, "deltaEstep.npy"))
            can_load = deltaEmin <= self.deltaEmin and deltaEmax >= self.deltaEmax and deltaEstep <= self.deltaEstep
        else:
            os.mkdir(self.dirName, 744)
        if can_load:
            self.deltaEstep = deltaEstep
            self.deltaEmax = deltaEmax
            self.deltaEmin = deltaEmin
            self.deltaEvals = np.arange(deltaEmin, deltaEmax, deltaEstep)
            self.vals = np.load(os.path.join(self.dirName, "vals.npy"))
        else:
            self.deltaEvals = np.arange(self.deltaEmin, self.deltaEmax, self.deltaEstep)
            self.vals = self.rateFunc(self.deltaEvals, self.Ec, self.otherParam, self.T)
            self.saveVals()
        return True

    def saveVals(self):
        np.save(os.path.join(self.dirName, "deltaEmin"), self.deltaEmin)
        np.save(os.path.join(self.dirName, "deltaEmax"), self.deltaEmax)
        np.save(os.path.join(self.dirName, "deltaEstep"), self.deltaEstep)
        np.save(os.path.join(self.dirName, "vals"), self.vals)

    def set_approx(self):
        self.approx = interp1d(self.deltaEvals, self.vals, assume_sorted=True)

    def increase_high_limit(self):
        new_deltaEmax = self.deltaEmax + np.abs(self.deltaEmax)
        new_inputs = np.arange(self.deltaEmax, new_deltaEmax, self.deltaEstep)
        new_vals = self.rateFunc(new_inputs,self.Ec, self.otherParam, self.T)
        self.deltaEmax = new_deltaEmax
        self.deltaEvals = np.hstack((self.deltaEvals, new_inputs))
        self.vals = np.hstack((self.vals, new_vals))
        self.saveVals()
        self.set_approx()
        print("High limit increased, Emax= " + str(self.deltaEmax))

    def decrease_low_limit(self):
        new_deltaEmin = self.deltaEmin - np.abs(self.deltaEmin)
        new_inputs = np.arange(new_deltaEmin, self.deltaEmin, self.deltaEstep)
        new_vals = self.rateFunc(new_inputs,self.Ec, self.otherParam, self.T)
        self.deltaEmin = new_deltaEmin
        self.deltaEvals = np.hstack((new_inputs, self.deltaEvals))
        self.vals = np.hstack((new_vals, self.vals))
        self.saveVals()
        self.set_approx()
        print("Low limit dencreased, Emin= " + str(self.deltaEmin))

    def get_tunnling_rates(self, deltaE):
        while self.deltaEmin - np.min(deltaE) >= -EPS:
            self.decrease_low_limit()
        while self.deltaEmax - np.max(deltaE) <= EPS:
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

    def __init__(self, rows, columns, VL, VR, VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, temperature, fastRelaxation=False,
                 tauLeaping=False, modifyR=False, constQ=False):
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
        self.constQ = constQ
        self.rows = rows
        self.columns = columns
        self.VL = VL
        self.VR = VR
        self.VG = flattenToColumn(VG)
        self.Q = flattenToColumn(Q0)
        self.n = flattenToColumn(n0)
        self.CG = flattenToColumn(CG)
        self.RG = flattenToColumn(RG)
        self.Ch = Ch
        self.Cv = np.vstack((np.zeros((1,columns)),Cv)) if Cv.size else np.zeros((1,1))
        self.Rh = Rh
        self.Rv = Rv if Rv.size else np.zeros((0,0))
        self.R = np.hstack((Rh.flatten(), Rh.flatten(), Rv.flatten(), Rv.flatten()))
        self.temperature = temperature
        self.variableRightWorkLen = (self.columns + 1) * self.rows
        self.variableDownWorkLen = self.columns * (self.rows - 1)
        self.totalChargePassedRight = 0
        self.totalChargePassedLeft = 0
        self.no_tunneling_next_time = False
        self.createCapacitanceMatrix()
        self.setDiagonalizedJ()
        self.setConstWork()
        self.setConstMatrix()
        self.setConstNprimePart()
        self.saveCurrentMap = False
        self.Ih = np.zeros(self.Rh.shape)
        self.Iv = np.zeros(self.Rv.shape)
        self.use_modifyR=modifyR
        # for std calculations
        # self.I_left_sqr_avg = 0
        # self.I_right_sqr_avg = 0

        # for variable R calculation
        if self.use_modifyR:
            self._leftnExponent = np.ones((self.rows, self.columns+1))
            self._rightnExponent = np.ones((self.rows, self.columns + 1))
        # Memory initialization for tau leaping
        self.tauLeaping = tauLeaping
        if tauLeaping:
            self.totalAction = np.zeros((self.rows, self.columns))
            self.right = np.zeros((self.rows, self.columns + 1))
            self.left = np.zeros((self.rows, self.columns + 1))
            self.down = np.zeros((self.rows - 1, self.columns))
            self.up = np.zeros((self.rows - 1, self.columns))

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
        copy_array.R = np.copy(self.R)
        copy_array.temperature = self.temperature
        copy_array.variableRightWorkLen = self.variableRightWorkLen
        copy_array.variableDownWorkLen = self.variableDownWorkLen
        copy_array.totalChargePassedRight = self.totalChargePassedLeft
        copy_array.totalChargePassedLeft = self.totalChargePassedRight
        copy_array.no_tunneling_next_time = self.no_tunneling_next_time
        copy_array.invC = self.invC
        copy_array._JeigenVectors = self._JeigenVectors
        copy_array._JeigenValues = self._JeigenValues
        copy_array._JeigenVectorsInv = self._JeigenVectorsInv
        copy_array.timeStep = self.timeStep
        copy_array.default_dt = self.default_dt
        copy_array._constQnPart = self._constQnPart
        copy_array._matrixQnPart = self._matrixQnPart
        copy_array.commonHorz = self.commonHorz
        copy_array.commonVert = self.commonVert
        copy_array.additionalLeft = self.additionalLeft
        copy_array.leftConstWork = self.leftConstWork
        copy_array.rightConstWork = self.rightConstWork
        copy_array.vertConstWork = self.vertConstWork
        copy_array.variableWork = np.copy(self.variableWork)
        copy_array.horizontalMatrix = self.horizontalMatrix
        copy_array.verticalMatrix = self.verticalMatrix
        copy_array._left_part_n_prime = self._left_part_n_prime
        copy_array._right_part_n_prime = self._right_part_n_prime
        copy_array.saveCurrentMap = self.saveCurrentMap
        copy_array.Ih = np.copy(self.Ih)
        copy_array.Iv = np.copy(self.Iv)
        # copy_array.I_left_sqr_avg = self.I_left_sqr_avg
        # copy_array.I_right_sqr_avg = self.I_right_sqr_avg
        copy_array.tauLeaping = self.tauLeaping
        copy_array.use_modifyR = self.use_modifyR
        copy_array._Q_eigenbasis = np.copy(self._Q_eigenbasis)
        copy_array.constQ = self.constQ
        if self.use_modifyR:
            copy_array._leftnExponent = np.copy(self._leftnExponent)
            copy_array._rightnExponent = np.copy(self._rightnExponent)
        copy_array.rates = np.copy(self.rates)
        copy_array.cumRates = np.copy(self.cumRates)
        # tau leaping
        if copy_array.tauLeaping:
            copy_array.totalAction = np.copy(self.totalAction)
            copy_array.left = np.copy(self.left)
            copy_array.right = np.copy(self.right)
            copy_array.down = np.copy(self.down)
            copy_array.up = np.copy(self.up)
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
        # self.I_left_sqr_avg = 0
        # self.I_right_sqr_avg = 0
        if self.saveCurrentMap:
            self.Ih[:,:]= 0
            self.Iv[:,:] = 0
        return True

    def changeVext(self, newVL, newVR):
        self.VL = newVL
        self.VR = newVR
        self.setConstWork()
        return True

    def getOccupation(self):
        return np.copy(self.n).reshape((self.rows, self.columns))

    def getGroundCharge(self):
        return np.copy(self.Q).reshape((self.rows, self.columns))

    def setOccupation(self, n):
        self.n = flattenToColumn(n)

    def setGroundCharge(self, Q):
        self.Q = flattenToColumn(Q)

    def getCharge(self):
        return self.totalChargePassedRight, self.totalChargePassedLeft

    # def getIsqr(self):
    #     return self.I_right_sqr_avg, self.I_left_sqr_avg

    def getTimeStep(self):
        return self.timeStep

    def getCurrentMap(self):
        if self.saveCurrentMap:
            map = np.zeros((self.rows*2 - 1, self.columns+1))
            map[::2,:] = np.copy(self.Ih)
            if self.Iv.size:
                map[1::2,:-1] = np.copy(self.Iv)
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
                -np.diagflat(n_diagonal, k=self.columns) - np.diagflat(n_diagonal, k=-self.columns)
        self.invC = np.linalg.inv(C_mat)
        return True

    def setConstNprimePart(self):
        self._left_part_n_prime = np.copy(self.Ch[:, :-1])
        self._left_part_n_prime[:, 1:] = 0
        self._right_part_n_prime = np.copy(self.Ch[:, 1:])
        self._right_part_n_prime[:, :-1] = 0
        self._right_part_n_prime = flattenToColumn(self._right_part_n_prime)
        self._left_part_n_prime = flattenToColumn(self._left_part_n_prime)
        return True

    def getNprime(self):
        return self.n + self._left_part_n_prime*self.VL + self._right_part_n_prime*self.VR

    def getNprimeForGivenN(self,n):
        return flattenToColumn(n) + self._left_part_n_prime*self.VL + self._right_part_n_prime*self.VR

    def getbVector(self):
        """
        For developement of Q by dQ/dt=JQ +b
        """
        res = -self.invC.dot(self.getNprime()) + self.VG
        return res / self.RG

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
        self.timeStep = -10/np.max(self._JeigenValues)
        self.default_dt = -0.1/np.min(self._JeigenValues)
        invMat = np.linalg.inv(self.invC + np.diagflat(1/self.CG))
        self._constQnPart = invMat.dot(self.VG)
        CGMat = np.repeat(1/flattenToRow(self.CG), invMat.shape[0], axis=0)
        self._matrixQnPart = invMat * CGMat - np.eye(invMat.shape[0])
        self._Q_eigenbasis = self._JeigenVectorsInv.dot(self.Q)
        return True

    def developeQ(self, dt):
        b = self._JeigenVectorsInv.dot(self.getbVector())
        exponent = np.exp(self._JeigenValues*dt)
        self._Q_eigenbasis = self._Q_eigenbasis*exponent +(b/self._JeigenValues)*(exponent - 1)
        self.Q = self._JeigenVectors.dot(self._Q_eigenbasis)
        return True

    def get_steady_Q_for_n(self):
        return self._constQnPart + self._matrixQnPart.dot(self.getNprime())

    def get_steady_Q_for_given_n(self, n):
        return self._constQnPart + self._matrixQnPart.dot(self.getNprimeForGivenN(n))

    def get_dist_from_steady(self, n, Q):
        return (flattenToColumn(Q) - self.get_steady_Q_for_given_n(n))**2

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
        self.leftConstWork = (0.5*self.commonHorz + self.additionalLeft).flatten()

        additionalRight = -self.additionalLeft
        self.rightConstWork = (0.5*self.commonHorz + additionalRight).flatten()

        self.vertConstWork = (0.5*self.commonVert).flatten()
        self.variableWork = np.zeros((4*self.rows*self.columns + 2*(self.rows - self.columns),))
        self.rates = np.zeros(self.variableWork.shape)
        self.cumRates = np.zeros(self.variableWork.shape)
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
        q = self.getNprime() + self.Q
        variableRightWork = (self.horizontalMatrix.dot(q)).flatten()
        variableDownWork = (self.verticalMatrix.dot(q)).flatten()
        variableRightWorkLen = variableRightWork.size
        variableDownWorkLen = variableDownWork.size
        self.variableWork[:variableRightWorkLen] = variableRightWork + self.rightConstWork
        self.variableWork[variableRightWorkLen:2*variableRightWorkLen] = -variableRightWork + self.leftConstWork
        self.variableWork[2 * variableRightWorkLen:2 * variableRightWorkLen + variableDownWorkLen] = variableDownWork + self.vertConstWork
        self.variableWork[2 * variableRightWorkLen + variableDownWorkLen:] = -variableDownWork + self.vertConstWork
        return self.variableWork

    def getRates(self):
        """
        Returns the tunnelling rate between neighboring dots
        :param T: temperature (in energy units = e^2[V])
        """
        if not self.constQ:
            self.getWork()
        work = np.copy(self.variableWork)
        if self.temperature == 0:
            work[work > -EPS] = 0
        else:
            exp_arg = work/self.temperature
            work[exp_arg > 20] = 0
            work[np.abs(exp_arg) < EPS] = -self.temperature
            rest = np.logical_and(EPS <= np.abs(exp_arg), np.abs(exp_arg)<= 20)
            work[rest] = work[rest]/(1-np.exp(exp_arg[rest]))
        if self.use_modifyR:
            self.modifyR()
        self.rates = -work / self.R
        return self.rates

    def getCurrentFromRates(self):
        current = np.sum(self.rates[:self.variableRightWorkLen]) -\
                      np.sum(self.rates[self.variableRightWorkLen:2*self.variableRightWorkLen])
        return current/(self.columns+1)

    def updateWorkForConstQ(self, fromDot, toDot):
        row1 = fromDot[0]*self.columns + fromDot[1] if 0 <= fromDot[1] < self.columns else None
        row2 = toDot[0]*self.columns + toDot[1] if 0 <= toDot[1] < self.columns else None
        variableRightWorkLen = (self.columns + 1)*self.rows
        variableDownWorkLen = self.columns*(self.rows - 1)
        if row1 is None:
            deltaRightWork = self.horizontalMatrix[:,row2]
            deltaDownWork = self.verticalMatrix[:,row2]
        elif row2 is None:
            deltaRightWork = - self.horizontalMatrix[:, row1]
            deltaDownWork = - self.verticalMatrix[:, row1]
        else:
            deltaRightWork = self.horizontalMatrix[:, row2] - self.horizontalMatrix[:, row1]
            deltaDownWork = self.verticalMatrix[:,row2] - self.verticalMatrix[:, row1]
        deltaLeftWork = -deltaRightWork

        deltaUpWork = -deltaDownWork
        self.variableWork[:variableRightWorkLen] += deltaRightWork
        self.variableWork[variableRightWorkLen:2 * variableRightWorkLen] += deltaLeftWork
        self.variableWork[2 * variableRightWorkLen:2 * variableRightWorkLen + variableDownWorkLen] += deltaDownWork
        self.variableWork[2 * variableRightWorkLen + variableDownWorkLen:] += deltaUpWork
        return self.variableWork

    def modifyR(self):
        nExponent = np.exp(-INV_DOS*self.n).reshape((self.rows, self.columns))
        self._rightnExponent[:,1:] = nExponent
        self._leftnExponent[:,:-1] = nExponent
        horzSize = self.Rh.size
        vertSize = self.Rv.size
        self.R[:horzSize] = (self.Rh*self._rightnExponent).flatten()
        self.R[horzSize:2*horzSize] = (self.Rh * self._leftnExponent).flatten()
        self.R[2*horzSize:2*horzSize+vertSize] = (self.Rv*nExponent[:-1,:]).flatten()
        self.R[2 * horzSize + vertSize:] = (self.Rv * nExponent[1:, :]).flatten()

    def getVoltages(self):
        return self.invC.dot(self.getNprime() + self.Q).reshape(self.rows, self.columns)

    def getVoltagesFromGround(self):
        return (self.VG - self.Q/self.CG).reshape(self.rows, self.columns)

    def getTimeInterval(self, randomNumber):
        """
        Calculates dt using Gillespie algorithm
        :return:
        """
        rates = self.getRates()
        np.cumsum(rates, out=self.cumRates)
        sum_rates = self.cumRates[-1]
        intervals = 0
        if sum_rates == 0:
            self.no_tunneling_next_time = True
            return self.default_dt, intervals
        dt = np.log(1 / randomNumber) / sum_rates
        if self.constQ:
            out_of_interval = False
        else:
            out_of_interval = dt > self.default_dt
        while out_of_interval:
            intervals += 1
            self.developeQ(self.default_dt)
            rates = self.getRates()
            np.cumsum(rates, out=self.cumRates)
            new_sum_rates = self.cumRates[-1]
            if new_sum_rates == 0:
                self.no_tunneling_next_time = True
                return self.default_dt, intervals
            dt = (dt - self.default_dt)*(new_sum_rates/sum_rates)
            sum_rates = new_sum_rates
            out_of_interval = dt > self.default_dt
            if np.max(self.get_dist_from_steady(self.n, self.Q)) < ALLOWED_ERR:
                break
        return dt, intervals

    def getLeapingTimeInterval(self):
        rates = self.getRates()
        if (rates <= EPS).all():
            self.no_tunneling_next_time = True
            return self.default_dt
        horzSize = self.rows * (self.columns + 1)
        vertSize = (self.rows - 1) * self.columns
        self.right[:,:] = rates[:horzSize].reshape(self.right.shape)
        self.left[:,:] = rates[horzSize:2*horzSize].reshape(self.left.shape)
        self.down[:,:] = rates[2*horzSize:2*horzSize + vertSize].reshape(self.up.shape)
        self.up[:,:] = rates[2*horzSize + vertSize:].reshape(self.down.shape)
        tunnelingTo = self.right[:,:-1] + self.left[:,1:]
        tunnelingTo[1:,:] += self.down
        tunnelingTo[:-1,:] += self.up
        tunnelingFrom = self.right[:,1:] + self.left[:,:-1]
        tunnelingFrom[1:,:] += self.up
        tunnelingFrom[:-1,:] += self.down
        changeAvg = np.abs(tunnelingTo - tunnelingFrom)
        changeVar = tunnelingTo + tunnelingFrom
        smallestChange = TAU_EPS*np.abs(self.getOccupation())
        smallestChange[smallestChange < 1] = 1
        tau = 0
        if (changeAvg > 0).any():
            tau = np.min(smallestChange[changeAvg > 0] / changeAvg[changeAvg > 0])
        if (changeVar > 0).any():
            tau2 = np.min(smallestChange[changeVar > 0]**2/changeVar[changeVar > 0])
            if tau > 0:
                tau = min(tau, tau2)
            else:
                tau = tau2
        if tau <= 0:
            tau = self.default_dt
            self.no_tunneling_next_time = True
        return tau

    def nextStep(self, dt, randomNumber):
        if not self.constQ:
            self.developeQ(dt)
        if self.no_tunneling_next_time:
            self.no_tunneling_next_time = False
            return True
        if (self.rates == 0).all():
            return True
        actionInd = np.searchsorted(self.cumRates, randomNumber*self.cumRates[-1])
        self.executeAction(actionInd)
        return True

    def nextLeapingSteps(self, dt):
        if not self.constQ:
            self.developeQ(dt)
        if self.no_tunneling_next_time:
            self.no_tunneling_next_time = False
            return True
        actionVec = np.random.poisson(lam=self.rates*dt, size=self.rates.size)
        self.executeMultipleActions(actionVec)

    def tunnel(self, fromDot, toDot, charge=1):
        if fromDot[1] == self.columns:
            self.totalChargePassedRight -= charge
        elif fromDot[1] == -1:
            self.totalChargePassedLeft -= charge
        else:
            self.n[fromDot[0]*self.columns + fromDot[1],0] -= charge
        if toDot[1] == self.columns:
            self.totalChargePassedRight += charge
        elif toDot[1] == -1:
            self.totalChargePassedLeft += charge
        else:
            self.n[toDot[0]*self.columns + toDot[1],0] += charge
        if self.constQ:
            self.updateWorkForConstQ(fromDot, toDot)
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

    def executeMultipleActions(self, actionVec, charge=1):
        horzSize = self.rows * (self.columns + 1)
        vertSize = (self.rows - 1) * self.columns
        self.right[:,:] = actionVec[:horzSize].reshape(self.right.shape)*charge
        self.left[:,:] = actionVec[horzSize:2 * horzSize].reshape(self.left.shape)*charge
        self.down[:,:] = actionVec[2 * horzSize:2 * horzSize + vertSize].reshape(self.down.shape)*charge
        self.up[:,:] = actionVec[2 * horzSize + vertSize:].reshape(self.up.shape)*charge
        total_passed_right = np.sum(self.right[:,-1] - self.left[:,-1])
        total_passed_left = np.sum(self.left[:,0] - self.right[:,0])
        self.totalChargePassedRight += total_passed_right
        self.totalChargePassedLeft += total_passed_left
        self.totalAction[:,:] = 0
        self.totalAction += self.left[:,1:] - self.left[:,:-1] + self.right[:,:-1] - self.right[:,1:]
        self.totalAction[1:,:] += self.down - self.up
        self.totalAction[:-1,:] += self.up - self.down
        self.n += flattenToColumn(self.totalAction)
        if self.saveCurrentMap:
            self.Ih += self.right - self.left
            self.Iv += self.down - self.up

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
        self.tunnel(fromDot, toDot, charge=charge)
        if self.constQ:
            self.updateWorkForConstQ(fromDot, toDot)
        return fromDot, toDot

class JJArray(DotArray):
    def __init__(self, rows, columns, VL, VR, VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv,
                 temperature, scGap, fastRelaxation=False, tauLeaping=False, modifyR=False):
        DotArray.__init__(self, rows, columns, VL, VR, VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv,
                          temperature, fastRelaxation=fastRelaxation, tauLeaping=tauLeaping, modifyR=modifyR)
        self.gap = scGap
        self.Ej = self.getEj()
        self.Ec = 1/np.mean(CG)
        self.qp_rate_calculator = TunnelingRateCalculator(-1, 1, 0.01, qp_tunneling, self.Ec, temperature, scGap,
                                                          "quasi_particles_rate")
        self.cp_rate_calculator = TunnelingRateCalculator(-1, 1, 0.01, cp_tunneling, self.Ec, temperature, self.Ej,
                                                          "cooper_pairs_rate")
        print("Rates were calculated")
        # dbg
        fig1 = self.qp_rate_calculator.plot_rate()
        fig2 = self.cp_rate_calculator.plot_rate()
        plt.show()
        plt.close(fig1)
        plt.close(fig2)

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
        leftConstWorkCp = 2*(self.commonHorz + self.additionalLeft).flatten()

        additionalRight = -self.additionalLeft
        rightConstWorkCp = 2*(self.commonHorz + additionalRight).flatten()

        vertConstWorkCp = (2*self.commonVert).flatten()

        self.constWorkCp = np.hstack((rightConstWorkCp, leftConstWorkCp, vertConstWorkCp, vertConstWorkCp))
        self.rates = np.zeros(2*self.variableWork.shape)

    def getWork(self):
        qp_work = super().getWork()
        cp_work = 2*self.variableWork + self.constWorkCp
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
            fromDot, toDot = super().executeAction(ind,charge=1)
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
    def __init__(self, index, VL0, VR0, Q0, n0, dotArray, constQ=False):
        self.dotArray = copy(dotArray)
        self.randomGenerator = np.random.RandomState()
        self.dotArray.setGroundCharge(Q0)
        self.dotArray.setOccupation(n0)
        self.VL = VL0
        self.VR = VR0
        self.index = index
        self.tauLeaping = dotArray.tauLeaping
        self.constQ = constQ
        if self.constQ:
            self.n = n0
            self.Q = Q0


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

    def executeLeapingStep(self, printState=False):
        dt = self.dotArray.getLeapingTimeInterval()
        self.dotArray.nextLeapingSteps(dt)
        if printState:
            self.printState()
        return dt

    def executeStep(self, printState=False):
        r = self.randomGenerator.rand(2)
        # for debug
        # r = next(self.randomGen)
        dt, intervals = self.dotArray.getTimeInterval(r[0])
        self.dotArray.nextStep(dt,r[1])
        if printState:
            self.printState()
        return dt + intervals*self.dotArray.default_dt

    def calcCurrent(self,print_stats=False, fullOutput=False, currentMap=False):
        if currentMap:
            self.dotArray.currentMapOn()
        if self.constQ:
            self.find_next_QG_using_gradient_descent()
            self.dotArray.setOccupation(self.n)
            self.dotArray.setGroundCharge(self.Q)
            self.dotArray.getWork()
        self.dotArray.resetCharge()
        final_t = self.dotArray.timeStep
        curr_t = 0
        steps = 0
        I_avg = 0
        I_var = 0
        if fullOutput:
            n_avg = np.zeros((self.dotArray.getRows(), self.dotArray.getColumns()))
            n_var = np.zeros(n_avg.shape)
            Q_avg = np.zeros((self.dotArray.getRows(), self.dotArray.getColumns()))
            Q_var = np.zeros(Q_avg.shape)
        curr_n = self.dotArray.getOccupation()
        curr_Q = self.dotArray.getGroundCharge()
        # err = 2*ALLOWED_ERR
        # errors = []
        while curr_t < final_t or steps < MIN_STEPS:
            if self.tauLeaping:
                dt = self.executeLeapingStep(printState=print_stats)
            else:
                dt = self.executeStep(printState=print_stats)
            steps += 1
            curr_I = self.dotArray.getCurrentFromRates()
            I_avg, I_var = self.update_statistics(curr_I, I_avg, I_var,
                                                  curr_t, dt)
            if fullOutput:
                n_avg, n_var = self.update_statistics(curr_n, n_avg, n_var, curr_t, dt)
                Q_avg, Q_var = self.update_statistics(curr_Q, Q_avg, Q_var, curr_t, dt)
                curr_n = self.dotArray.getOccupation()
                curr_Q = self.dotArray.getGroundCharge()
            curr_t += dt
            # if steps%MIN_STEPS == 0:
                # rightCurrent, leftCurrent = self.dotArray.getCharge()
                # rightCurrentSqr, leftCurrentSqr = self.dotArray.getIsqr()
                # rightCurrentErr = np.sqrt((1/(steps-1)) * (rightCurrentSqr/curr_t - rightCurrent**2/curr_t**2))
                # leftCurrentErr = np.sqrt((1/(steps-1)) * (leftCurrentSqr/curr_t - leftCurrent**2/curr_t**2))
                # err = max(rightCurrentErr, leftCurrentErr)
                # err = self.get_err(I_var, steps, curr_t)
                # errors.append(err)
                # if err < I_avg*0.01:
                #     break
        err = self.get_err(I_var, steps, curr_t)
        res = (I_avg, err)
        if fullOutput:
            res = res + (n_avg, Q_avg, self.get_err(n_var, steps, curr_t), self.get_err(Q_var, steps, curr_t))
        if currentMap:
            res = res + (self.dotArray.getCurrentMap()/curr_t,)
            self.dotArray.currentMapOff()
        return res

    def calcAverageNForGivenQ(self, Q, n0, calcVoltages=False):
        curr_t = 0
        steps = 0
        n_avg = np.zeros((self.dotArray.getRows(), self.dotArray.getColumns()))
        n_var = np.zeros(n_avg.shape)
        self.dotArray.setGroundCharge(Q)
        self.dotArray.setOccupation(n0)
        self.dotArray.getWork()
        curr_n = n0
        err = 2*ALLOWED_ERR
        while err > ALLOWED_ERR:
            if self.tauLeaping:
                dt = self.executeLeapingStep()
            else:
                dt = self.executeStep()
            steps += 1
            n_avg, n_var = self.update_statistics(curr_n, n_avg, n_var, curr_t, dt)
            curr_n = self.dotArray.getOccupation()
            curr_t += dt
            if steps%MIN_STEPS == 0:
                err = np.max(self.get_err(n_var, steps, curr_t))
        res = (n_avg, self.get_err(n_var, steps, curr_t))
        if calcVoltages:
            self.dotArray.setOccupation(n_avg)
            v = self.dotArray.getVoltages()
            v_from_ground = self.dotArray.getVoltagesFromGround()
            res = res + (v, v_from_ground)
        return res

    def calc_lyaponuv_grad(self, Q):
        n_avg, _ = self.calcAverageNForGivenQ(Q, self.n)
        self.dotArray.setOccupation(n_avg)
        self.dotArray.setGroundCharge(Q)
        self.n = np.round(n_avg)
        return Q.flatten() - self.dotArray.get_steady_Q_for_n().flatten()

    def find_next_QG_using_gradient_descent(self):
        flag = False
        rep = 0
        lr = INI_LR
        while (not flag) and rep < GRAD_REP:
            res, flag = simple_gadient_descent(self.calc_lyaponuv_grad, self.Q, lr=lr, plot_lc=False)
            rep += 1
            lr = lr/10
        if flag: #if gradient descent did not converge skipping the point
            self.Q = res

    def plotAverageVoltages(self, Q0, dQ, n0, points_num):
        size = self.dotArray.rows*self.dotArray.columns
        vs = np.zeros((points_num,size))
        vgs = np.zeros((points_num,size))
        Qmagnitude = np.zeros((points_num,))
        Q=Q0
        for point in range(points_num):
            _,_,v,vg = self.calcAverageNForGivenQ(Q, n0, calcVoltages=True)
            vs[point,:] = v
            vgs[point,:] = vg
            Qmagnitude[point] = np.sqrt(np.sum(Q**2))
            Q += dQ
        for dot in range(size):
            plt.figure()
            plt.plot(Qmagnitude, vs[:,dot], 'b',
                     Qmagnitude, vgs[:,dot], 'r')
        plt.show()

    def getToSteadyState(self):
        curr_t = 0
        steps = 0
        n_avg = np.zeros(
            (self.dotArray.getRows(), self.dotArray.getColumns()))
        n_var = np.zeros(n_avg.shape)
        Q_avg = np.zeros(n_avg.shape)
        Q_var =  np.zeros(n_avg.shape)
        curr_n = self.dotArray.getOccupation()
        curr_Q = self.dotArray.getGroundCharge()
        err = ALLOWED_ERR*2
        # plot = self.dotArray.VL - self.dotArray.VR > 4
        # if plot:
        #     Qs = []
        #     Qn = []
        #     ts = []
        #     t=0
        while err > ALLOWED_ERR:
            if self.tauLeaping:
                dt = self.executeLeapingStep()
            else:
                dt = self.executeStep()
            steps += 1
            n_avg, n_var = self.update_statistics(curr_n, n_avg, n_var,
                                                  curr_t, dt)
            Q_avg, Q_var = self.update_statistics(curr_Q, Q_avg, Q_var,
                                                  curr_t, dt)
            curr_n = self.dotArray.getOccupation()
            curr_Q = self.dotArray.getGroundCharge()
            curr_t += dt
            # if plot:
            #     t+=dt
            #     Qs.append(Q_avg)
            #     Qn.append(self.dotArray.get_steady_Q_for_given_n(n_avg).reshape(Qs[0].shape))
            #     ts.append(t)
            if steps % MIN_STEPS == 0:
                new_err = np.max(self.dotArray.get_dist_from_steady(n_avg, Q_avg))
                err = new_err
                n_avg = np.zeros(
                    (self.dotArray.getRows(), self.dotArray.getColumns()))
                n_var = np.zeros(n_avg.shape)
                Q_avg = np.zeros(n_avg.shape)
                Q_var = np.zeros(n_avg.shape)
                curr_t = 0
                # if plot and steps % (100 * MIN_STEPS):
                #     Qs = np.array(Qs)
                #     Qn = np.array(Qn)
                #     for i in range(Qs.shape[2]):
                #         plt.plot(ts, Qs[:,:,i],'.')
                #         plt.plot(ts, Qn[:, :, i], '*')
                #     print(err)
                #     plt.show()
                #     Qs = []
                #     Qn = []
                #     ts = []

        return True

    def checkSteadyState(self, rep):
        I = []
        ns = []
        Qs = []
        t = self.dotArray.getTimeStep()/10
        for r in range(rep):
            # running once to get to steady state
            self.calcCurrent(t)
            # now we are in steady state calculate current
            current, currentErr, n, Q, nErr, QErr =\
                self.calcCurrent(t, fullOutput=True)
            I.append(current)
            ns.append(n)
            Qs.append(Q)
        ns = np.array(ns)
        Qs = np.array(Qs)
        x = t*np.arange(rep)
        print("Mean: " + str(np.average(I))+ " Var right: " + str(np.var(I)))
        print("Tail mean : " + str(np.average(I[len(I)//4:])) + " Tail var right: " + str(np.var(I[len(I)//4:])))
        plt.figure()
        plt.plot(x,I)
        plt.figure()
        plt.plot(x, ns.reshape((ns.shape[0], ns.shape[1]*ns.shape[2])))
        plt.figure()
        plt.plot(x, Qs.reshape((Qs.shape[0], Qs.shape[1]*Qs.shape[2])))
        plt.show()

    def saveState(self, I, IErr, n=None, Q=None, nErr=None, QErr=None, Imaps=None,
                  fullOutput=False, currentMap=False, basePath=''):
        baseName = basePath + "_temp_" + str(self.index)
        if fullOutput:
            np.save(baseName + "_ns", np.array(n))
            np.save(baseName + "_Qs", np.array(Q))
            np.save(baseName + "_nsErr", np.array(nErr))
            np.save(baseName + "_QsErr", np.array(QErr))
        np.save(baseName + "_n", self.dotArray.getOccupation())
        np.save(baseName + "_Q", self.dotArray.getGroundCharge())
        if currentMap:
            np.save(baseName + "_current_map", np.array(Imaps))
        np.save(baseName + "_I", np.array(I))
        np.save(baseName + "_IErr", np.array(IErr))

    def loadState(self,  fullOutput=False, currentMap=False, basePath=''):
        baseName = basePath + "_temp_" + str(self.index)
        I = np.load(baseName + "_I.npy")
        loadLen = len(I)
        IErr = np.load(baseName + "_IErr.npy")
        if len(IErr) > loadLen:
            IErr = IErr[:loadLen]
        n = np.load(baseName + "_n.npy")
        Q = np.load(baseName + "_Q.npy")
        res = (I,IErr,n,Q)
        if fullOutput:
            ns = np.load(baseName + "_ns.npy")
            if len(ns) > loadLen:
                ns = ns[:, loadLen, :, :]
            Qs = np.load(baseName + "_Qs.npy")
            if len(Qs) > loadLen:
                Qs = Qs[:, loadLen, :, :]
            nsErr = np.load(baseName + "_nsErr.npy")
            if len(nsErr) > loadLen:
                nsErr = nsErr[:, loadLen, :, :]
            QsErr = np.load(baseName + "_QsErr.npy")
            if len(QsErr) > loadLen:
                QsErr = QsErr[:, loadLen, :, :]
            res = res + (ns, Qs,nsErr, QsErr)
        if currentMap:
            Imaps = np.load(baseName + "_current_map.npy")
            if len(Imaps) > loadLen:
                Imaps = Imaps[:, loadLen, :, :]
            res = res + (Imaps,)
        return res

    def calcIV(self, Vmax, Vstep, vSym, fullOutput=False, print_stats=False,
               currentMap=False, basePath="", resume=False):
        I = []
        IErr = []
        ns = []
        Qs = []
        nsErr = []
        QsErr = []
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
        if resume:
            resumeParams = self.loadState(fullOutput=fullOutput, currentMap=currentMap, basePath=basePath)
            I = list(resumeParams[0])
            IErr = list(resumeParams[1])
            Vind = len(I)
            VL_vec = VL_vec[Vind:]
            VR_vec = VR_vec[Vind:]
            self.dotArray.setOccupation(resumeParams[2])
            self.dotArray.setGroundCharge(resumeParams[3])
            if fullOutput:
                ns = list(resumeParams[4])
                Qs = list(resumeParams[5])
                nsErr = list(resumeParams[6])
                QsErr = list(resumeParams[7])
            if currentMap:
                Imaps = list(resumeParams[8])
        for VL,VR in zip(VL_vec, VR_vec):
            self.dotArray.changeVext(VL, VR)
            # running once to get to steady state
            if not self.constQ:
                self.getToSteadyState()
            # now we are in steady state calculate current
            stepRes = self.calcCurrent(print_stats=print_stats, fullOutput=fullOutput, currentMap=currentMap)
            current = stepRes[0]
            currentErr = stepRes[1]
            if fullOutput:
                ns.append(stepRes[2])
                Qs.append(stepRes[3])
                nsErr.append(stepRes[4])
                QsErr.append(stepRes[5])
            if currentMap:
                Imaps.append(stepRes[-1])
            I.append(current)
            IErr.append(currentErr)
            if self.index == 0:
                print(VL - VR, end=',')
            self.saveState(I, IErr, ns, Qs, nsErr, QsErr, Imaps, fullOutput=fullOutput,
                           currentMap=currentMap, basePath=basePath)
        res = (np.array(I), np.array(IErr), VL_res - VR_res)
        if fullOutput:
            res = res + (np.array(ns), np.array(Qs), np.array(nsErr), np.array(QsErr))
        if currentMap:
            res = res + (np.array(Imaps),)
        return res

    def printState(self):
        self.dotArray.printState()

    def update_statistics(self, value, avg, n_var, total_time, time_step):
        """
        Updating the statistics of a measured value according to West's
        algorithm (as described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2973983/)
        :param value: Measured value this step
        :param avg: Average up to this step (result of last call to this function)
        :param n_var: non-normalized variance up to this step (result of last call to this function)
        :param total_time: Total simulation time not including last step
        :param new_time: The time of the last step
        :return: (avg, n_var) Upadted values.
        to get the actual variance use Var = (n*n_var)/((n-1)*total_time)
        """
        new_time = total_time + time_step
        dist_from_avg = value - avg
        local_std = dist_from_avg*time_step / new_time
        new_n_var = n_var + dist_from_avg*total_time*local_std
        new_avg = avg + local_std
        return new_avg, new_n_var

    def get_err(self, n_var, steps, time):
        """
        Returns the standard error of mean of a variable based on the
        non-normalized variance (see method "update_statistics")
        :param n_var: non-normalized variance
        :param steps: number of steps for data point
        :param time: total running time for data point
        :return: standart deviation
        """
        return np.sqrt((1/((steps-1)*time))*n_var)

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
        self.dotArray.setGroundCharge(np.copy(Q))
        states.append(np.copy(self.n0).flatten())
        tup_n = tuple(self.n0.flatten())
        states_dict[tup_n] = 0
        left_rates_diff = []
        right_rates_diff = []
        current_state_ind = 0
        next_state_ind = 1
        edges_line = [0]
        while current_state_ind < len(states):
            self.dotArray.setOccupation(self.reshape_to_array(np.copy(states[current_state_ind])))
            rates = self.dotArray.getRates()
            for ind,rate in enumerate(rates):
                if rate > 0: # add new edge
                    fromDot, toDot = self.dotArray.executeAction(ind,1)
                    n = self.dotArray.getOccupation().flatten()
                    tup_n = tuple(n)
                    if tup_n not in states_dict:
                        states_dict[tup_n] = next_state_ind
                        states.append(n)
                        next_state_ind += 1
                        edges_line.append(rate)
                    else:
                        edges_line[states_dict[tup_n]] += rate
                    self.dotArray.tunnel(toDot, fromDot,1)
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
        n_copy = self.dotArray.getOccupation()
        Q_copy = self.dotArray.getGroundCharge()
        n = self.get_average_state(Q)
        self.dotArray.setOccupation(n)
        self.dotArray.setOccupation(Q)
        v = self.dotArray.getVoltages()
        self.dotArray.setOccupation(n_copy)
        self.dotArray.setGroundCharge(Q_copy)
        return v

    def get_voltages_from_ground(self, Q):
        Q_copy = self.dotArray.getGroundCharge()
        self.dotArray.setGroundCharge(Q)
        v = self.dotArray.getVoltagesFromGround()
        self.dotArray.setGroundCharge(Q_copy)
        return v

    def calc_lyaponuv(self, Q):
        n_copy = self.dotArray.getOccupation()
        Q_copy = self.dotArray.getGroundCharge()
        n = self.get_average_state(Q)
        self.dotArray.setOccupation(n)
        self.dotArray.setGroundCharge(Q)
        diff_from_equi = np.sum(flattenToColumn(Q) - self.dotArray.get_steady_Q_for_n())
        self.dotArray.setOccupation(n_copy)
        self.dotArray.setGroundCharge(Q_copy)
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
        self.dotArray.setOccupation(n)
        self.dotArray.setGroundCharge(Q)
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
        # self.find_next_QG_using_lyaponuv(basePath)
        # print(self.QG)
        # self.plot_average_voltages(self.QG-Q_SHIFT, self.QG+Q_SHIFT)
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

    def saveState(self, I, n=None, Q=None, fullOutput=False, basePath=''):
        baseName = basePath + "_temp_" + str(self.index)
        if fullOutput:
            np.save(baseName + "_ns", np.array(n))
            np.save(baseName + "_Qs", np.array(Q))
        np.save(baseName + "_n", self.n0)
        np.save(baseName + "_Q", self.QG)
        np.save(baseName + "_I", np.array(I))


    def loadState(self, fullOutput=False, basePath=''):
        baseName = basePath + "_temp_" + str(self.index)
        I = np.load(baseName + "_I.npy")
        n = np.load(baseName + "_n.npy")
        Q = np.load(baseName + "_Q.npy")
        res = (I, n, Q)
        if fullOutput:
            ns = np.load(baseName + "_ns.npy")
            Qs = np.load(baseName + "_Qs.npy")
            res = res + (ns, Qs)
        return res

    def calcIV(self, Vmax, Vstep, vSym, fullOutput=False, print_stats=False, currentMap=False,
               basePath="", resume=False):
        # TODO: add error calculation
        I = []
        ns = []
        Qs = []
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
            Vind = len(I)
            VL_vec = VL_vec[Vind:]
            VR_vec = VR_vec[Vind:]
            self.n0 = resumeParams[1]
            self.QG = resumeParams[2]
            if fullOutput:
                ns = list(resumeParams[3])
                Qs = list(resumeParams[4])
        for VL,VR in zip(VL_vec,VR_vec):
            if self.index == 0:
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
            self.saveState(I, ns, Qs, fullOutput=fullOutput, basePath=basePath)
        result = (np.array(I), np.zeros((len(I),)), VL_res - VR_res)
        if fullOutput:
            result = result + (ns, Qs, np.zeros((len(ns),)), np.zeros((len(Qs),)))
        return result

def runSingleSimulation(index, VL0, VR0, vSym, Q0, n0,Vmax, Vstep, dotArray,
                        fullOutput=False, printState=False, useGraph=False, currentMap=False,
                        basePath="", resume=False, constQ=False):
    if useGraph:
        simulator = GraphSimulator(index, VL0, VR0, Q0, n0, dotArray)
    else:
        simulator = Simulator(index, VL0, VR0, Q0, n0, dotArray, constQ)
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
                      gap=0, leaping=False, modifyR=False, plotVoltages=False,
                      constQ=False):
    basePath = os.path.join(savePath, fileName)
    if useGraph:
        fastRelaxation = True
        repeats = 1
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
                                 temperature, gap, fastRelaxation=fastRelaxation, tauLeaping=leaping, modifyR=modifyR)
        print("Superconducting prototype array was created")
    else:
        prototypeArray = DotArray(rows, columns, VL0, VR0, np.array(VG0), np.array(Q0), np.array(n0), np.array(CG),
                                  np.array(RG), np.array(Ch), np.array(Cv), np.array(Rh), np.array(Rv),
                                  temperature, fastRelaxation=fastRelaxation, tauLeaping=leaping, modifyR=modifyR,
                                  constQ=constQ)
        print("Normal prototype array was created")
    if checkSteadyState:
        print("Running steady state convergance check")
        simulator = Simulator(0, VL0, VR0, Q0, n0, prototypeArray)
        simulator.checkSteadyState(repeats)
        simulator.dotArray.changeVext(VL0+Vstep, VR0)
        simulator.checkSteadyState(repeats)
        exit(0)

    if plotVoltages:
        print("Plotting voltages")
        points_num = 100
        Q0 = np.array([[1,1,0.1]])
        dQ = np.array([[0,0,0.01]])
        VL0 = 1
        simulator = Simulator(0, VL0, VR0, Q0, n0, prototypeArray)
        simulator.plotAverageVoltages(Q0,dQ,n0,points_num)
        exit(0)
    Is = []
    IsErr = []
    ns = []
    Qs = []
    nsErr = []
    QsErr = []
    Imaps = []
    if not dbg:
        print("Starting parallel run")
        pool = Pool(processes=repeats)
        results = []
        for repeat in range(repeats):
            res = pool.apply_async(runSingleSimulation,
                                    (repeat, VL0, VR0, vSym, Q0, n0, Vmax, Vstep, prototypeArray, fullOutput,
                                     printState, useGraph, currentMap,basePath, resume, constQ))
            results.append(res)
        for res in results:
            result = res.get()
            I = result[0]
            IErr = result[1]
            currentMapInd = 3
            if fullOutput:
                n = result[3]
                Q = result[4]
                nErr = result[5]
                QErr = result[6]
                ns.append(n)
                Qs.append(Q)
                nsErr.append(nErr)
                QsErr.append(QErr)
                currentMapInd = 7
            if currentMap:
                Imaps.append(result[currentMapInd])
            Is.append(I)
            IsErr.append(IErr)
        V = result[2]
        params = result[-1]

    else: #dbg
        print("Starting serial run")
        for repeat in range(repeats):
            result = runSingleSimulation(repeat, VL0, VR0, vSym, Q0, n0, Vmax, Vstep, prototypeArray, fullOutput,
                                         printState, useGraph, currentMap,basePath, resume, constQ)
            I = result[0]
            IErr = result[1]
            currentMapInd = 3
            if fullOutput:
                n = result[3]
                Q = result[4]
                nErr = result[5]
                QErr = result[6]
                ns.append(n)
                Qs.append(Q)
                nsErr.append(nErr)
                QsErr.append(QErr)
                currentMapInd = 7
            if currentMap:
                Imaps.append(result[currentMapInd])
            Is.append(I)
            IsErr.append(IErr)
        V = result[2]
        params = result[-1]
        # dbg
    print("Saving results")
    avgI = np.mean(np.array(Is), axis=0)
    if repeats < 10:
        avgIErr = np.sqrt(np.sum(np.array(IsErr)**2, axis=0)/len(IsErr))
    else:
        avgIErr = np.std(np.array(Is), axis=0)/np.sqrt(len(Is))
    if fullOutput:
        avgN = np.mean(np.array(ns), axis=0)
        avgQ = np.mean(np.array(Qs), axis=0)
        avgNErr = np.sqrt(np.sum(np.array(nsErr) ** 2, axis=0)) / len(nsErr)
        avgQErr = np.sqrt(np.sum(np.array(QsErr) ** 2, axis=0)) / len(QsErr)
        np.save(basePath + "_n", avgN)
        np.save(basePath + "_Q", avgQ)
        np.save(basePath + "_nErr", avgNErr)
        np.save(basePath + "_QErr", avgQErr)
        np.save(basePath + "_full_I", np.array(Is))
        np.save(basePath + "_full_IErr", np.array(IsErr))
    fig = plt.figure()
    IplusErr = avgI + avgIErr
    IminusErr = avgI - avgIErr
    plt.plot(V[:V.size//2], avgI[:V.size//2], 'b-',V[:V.size//2], IplusErr[:V.size//2], 'b--',
             V[:V.size//2], IminusErr[:V.size//2], 'b--',
             V[V.size//2:], avgI[V.size//2:], 'r-',V[V.size//2:], IplusErr[V.size//2:], 'r--',
             V[V.size//2:], IminusErr[V.size//2:], 'r--')
    plt.savefig(basePath + "_IV.png")
    np.save(basePath + "_I", avgI)
    np.save(basePath + "_IErr", avgIErr)
    np.save(basePath + "_V", V)
    plt.close(fig)
    if currentMap:
        avgImaps = np.mean(np.array(Imaps),axis=0)
        np.save(basePath + "_Imap", avgImaps)
    for index in range(repeats):
        removeState(index, fullOutput=fullOutput, basePath=basePath, currentMap=currentMap, graph=use_graph)
    removeRandomParams(basePath)
    return params

def removeState(index, fullOutput=False, basePath='', currentMap=False, graph=False):
    baseName = basePath + "_temp_" + str(index)
    os.remove(baseName + "_I.npy")
    if not graph:
        os.remove(baseName + "_IErr.npy")
    os.remove(baseName + "_n.npy")
    os.remove(baseName + "_Q.npy")
    if fullOutput:
        os.remove(baseName + "_ns.npy")
        os.remove(baseName + "_Qs.npy")
        if not graph:
            os.remove(baseName + "_nsErr.npy")
            os.remove(baseName + "_QsErr.npy")
    if currentMap:
        os.remove(baseName + "_current_map.npy")
    return True

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
    parser.add_option("--tau-leaping", dest="leaping", help="use tau leaping approximation [Default:%default]",
                      default=False, action='store_true')
    parser.add_option("--const-q", dest="constQ",
                      help="calc current using const Q method [Default:%default]",
                      default=False, action='store_true')
    parser.add_option("--variable-ef", dest="modifyR", help="if true Fermi energy level for each island will be changed "
                                                            "according to constant density of states assumption, else"
                                                            " it will be assumed constant (infinite density of states"
                                                            " [Default:%default]",
                      default=False, action='store_true')
    parser.add_option("-o", "--output-folder", dest="output_folder",
                      help="Output folder [default: current folder]",
                      default='.')
    parser.add_option("-i", "--load-from-file", dest="params_path",
                      help="If a parameters file is given all the array parameters would be loaded from that"
                           " file. The file should be in the same format as the resulted parameter file for a run."
                           " Ignores all other array related parameters, so the array owuld be exactly as specified in"
                           " the given file. [default: '']",
                      default='')


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

def load_params_from_file(file_path):
    runningParams = dict()
    arrayParams = dict()
    with open(file_path, 'r') as f:
        start = False
        for line in f:
            if ' = ' in line:
                splitted = line.split(' = ')
                key = splitted[0]
                try:
                    value = literal_eval(splitted[1].rstrip('\n'))
                except Exception:
                    value = splitted[1].rstrip('\n')
                runningParams[key] = value
            elif ': ' in line:
                if start:
                    key = splitted[0]
                    try:
                        splitted[1] = splitted[1].rstrip('\n')
                        splitted[1] = re.sub('\[\s+','[',splitted[1])
                        splitted[1] = re.sub('\s+\]', ']', splitted[1])
                        splitted[1] = re.sub('\s+', ',', splitted[1])
                        value = literal_eval(splitted[1])
                    except Exception:
                        value = splitted[1].rstrip('\n')
                    arrayParams[key] = value
                start = True
                splitted = line.split(': ')
            elif start:
                splitted[1] = splitted[1].replace('\n',' ') + line
        key = splitted[0]
        try:
            splitted[1] = splitted[1].rstrip('\n')
            splitted[1] = re.sub('\[\s+', '[', splitted[1])
            splitted[1] = re.sub('\s+\]', ']', splitted[1])
            splitted[1] = re.sub('\s+', ',', splitted[1])
            value = literal_eval(splitted[1])
        except Exception:
            value = splitted[1].rstrip('\n')
        arrayParams[key] = value
    return runningParams, arrayParams

if __name__ == "__main__":
    # Initializing Running Parameters
    options, args = getOptions()
    params_file = options.params_path
    if params_file:
        runningParams, arrayParams = load_params_from_file(params_file)
        rows = runningParams['M']
        columns = runningParams['N']
    else:
        rows = options.M
        columns = options.N
    VR0 = options.VR
    VL0 = VR0 + options.Vmin
    dist = options.dist
    T = options.T
    gap = options.gap
    sc = options.sc
    leaping = options.leaping
    Q0 = create_random_array(rows, columns, options.Q0_avg, options.Q0_std, dist,
                             False)
    n0 = create_random_array(rows, columns, options.n0_avg, options.n0_std, dist, False)
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
    modifyR = options.modifyR
    constQ = options.constQ
    if params_file:
        VG = arrayParams['VG']
        CG = arrayParams['CG']
        RG = arrayParams['RG']
        Ch = arrayParams['Ch']
        Cv = arrayParams['Cv']
        Cv = Cv[1:]
        Rh = arrayParams['Rh']
        Rv = arrayParams['Rv']
    else:
        VG = create_random_array(rows, columns, options.VG_avg, options.VG_std, dist,
                                 False)
        CG = create_random_array(rows, columns, options.CG_avg, options.CG_std, dist, True)
        RG = create_random_array(rows, columns, options.RG_avg, options.RG_std, dist, True)
        if options.custom_ch.replace('\"', '') and options.custom_cv.replace('\"', ''):
            Ch = literal_eval(options.custom_ch.replace('\"', ''))
            Cv = literal_eval(options.custom_cv.replace('\"', ''))
        else:
            Ch = create_random_array(rows, columns + 1, options.C_avg, options.C_std, dist,
                                     True)
            Cv = create_random_array(rows - 1, columns, options.C_avg, options.C_std, dist,
                                     True)
        if options.custom_rh.replace('\"', '') and options.custom_rv.replace('\"', ''):
            Rh = literal_eval(options.custom_rh.replace('\"', ''))
            Rv = literal_eval(options.custom_rv.replace('\"', ''))
        else:
            Rh = create_random_array(rows, columns + 1, options.R_avg, options.R_std, dist,
                                     True)
            Rv = create_random_array(rows - 1, columns, options.R_avg, options.R_std, dist,
                                     True)

    # Running Simulation
    if not os.path.exists(savePath):
        os.mkdir(savePath)
    elif not os.path.isdir(savePath):
        print("the given path exists but is a file")
        exit(0)
    if dbg:
        import cProfile, pstats, io
        from pstats import SortKey
        pr = cProfile.Profile()
        pr.enable()

    array_params = runFullSimulation(VL0, VR0, vSym, VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, rows,  columns,
                                     Vmax, Vstep, temperature=T, repeats=repeats, savePath=savePath, fileName=fileName,
                                     fullOutput=fullOutput, printState=False, useGraph=use_graph,
                                     fastRelaxation=fast_relaxation, currentMap=current_map,
                                     dbg=dbg, plotCurrentMaps=plot_current_map, resume=resume,
                                     checkSteadyState=False, superconducting=sc, gap=gap, leaping=leaping,
                                     modifyR=modifyR, constQ=constQ)
    saveParameters(savePath, fileName, options, array_params)

    if dbg:
        pr.disable()
        s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())

    exit(0)
