__author__ = 'shahar'

import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import numpy as np
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
import matplotlib
# matplotlib.use("Agg")
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import matplotlib.animation as animation
# from mayavi import mlab
from optparse import OptionParser
import re
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
MIN_STEPS = 100
STEADY_STATE_VAR = 1e-4
ALLOWED_ERR = 1e-4
STEADY_STATE_REP = 100
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
    exp_arg = x/T
    if exp_arg > 20:
        return 0
    else:
        return 1/(1 + np.exp(x/T))

def qp_density_of_states(x,energy_gap):
    arg = x**2-energy_gap**2
    if arg < EPS:
        arg = EPS
    return np.abs(x) / np.sqrt(arg)

def cp_tunneling(x, Ec, Ej, T):
    return (np.pi*Ej)**2 * high_impadance_p(x, Ec, T, 2)

def qp_integrand(x2, x1, deltaE, Ec, gap, T):
    return qp_density_of_states(x1, gap)*qp_density_of_states(x2, gap)*fermi_dirac_dist(-x1, T)\
           *fermi_dirac_dist(x2, T)*high_impadance_p(x2-x1 - deltaE, Ec, T, 1)



def qp_tunneling_single(deltaE, Ec, gap, T):
    part1 = dblquad(qp_integrand, -np.infty, -gap-EPS, lambda x: -np.infty,
                    lambda x: -gap-EPS, args=(deltaE,Ec, gap, T))[0]
    part2 = dblquad(qp_integrand, -np.infty, -gap-EPS, lambda x: gap,
                    lambda x: np.infty, args=(deltaE,Ec, gap, T))[0]
    part3 = dblquad(qp_integrand, gap+EPS, np.infty, lambda x: -np.infty,
                    lambda x: -gap-EPS, args=(deltaE,Ec, gap, T))[0]
    part4 = dblquad(qp_integrand, gap+EPS, np.infty, lambda x: gap,
                    lambda x: np.infty+EPS, args=(deltaE,Ec, gap, T))[0]
    return part1 + part2 + part3 + part4

def qp_tunneling(deltaE, Ec, gap, T):
        res = np.zeros(deltaE.shape)
        for ind,val in enumerate(deltaE):
            res[ind] = np.exp(val/T)*qp_tunneling_single(val, Ec, gap, T)
            # res[ind] =  qp_tunneling_single(val, Ec, gap, T)
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
        print("Increasing high limit")
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
        print("Decreasing low limit")
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
        while self.deltaEmin - np.min(deltaE) >= self.deltaEstep:
            self.decrease_low_limit()
        while self.deltaEmax - np.max(deltaE) <= -self.deltaEstep:
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

    def __init__(self, rows, columns, VL, VR, VU, VD, VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, temperature, fastRelaxation=False,
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
        self.rows = rows
        self.columns = columns
        self.VL = VL
        self.VR = VR
        self.VU = VU
        self.VD = VD
        self.VG = flattenToColumn(VG)
        self.n = flattenToColumn(n0)
        self.CG = flattenToColumn(CG)
        self.RG = flattenToColumn(RG)
        self.Ch = Ch
        self.Cv = Cv
        self.Rh = Rh
        self.Rv = Rv
        self.R = np.hstack((Rh.flatten(), Rh.flatten(), Rv.flatten(), Rv.flatten()))
        self.temperature = temperature
        self.rightWorkLen = (self.columns + 1) * self.rows
        self.downWorkLen = self.columns * (self.rows + 1)
        self.work = np.zeros((self.rightWorkLen*2 + self.downWorkLen*2,))
        self.no_tunneling_next_time = False
        self.createConversionMatrices()
        self.setInitialVoltages(Q0)
        self.setEquilibriumVoltages()
        self.setChangeMatrices()
        self.use_modifyR=modifyR
        self.rates = np.zeros((self.rightWorkLen*2 + self.downWorkLen*2,))
        self.cumRates = np.zeros(self.rates.shape)
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
            self.down = np.zeros((self.rows + 1, self.columns))
            self.up = np.zeros((self.rows + 1, self.columns))

    def __copy__(self):
        copy_array = object.__new__(type(self))
        copy_array.rows = self.rows
        copy_array.columns = self.columns
        copy_array.VL = self.VL
        copy_array.VR = self.VR
        copy_array.VU = self.VU
        copy_array.VD = self.VD
        copy_array.VG = self.VG
        copy_array.voltagesDiag= np.copy(self.voltagesDiag)
        copy_array.voltages = np.copy(self.voltages)
        copy_array.n = np.copy(self.n)
        copy_array.CG = self.CG
        copy_array.RG = self.RG
        copy_array.Ch = self.Ch
        copy_array.Cv = self.Cv
        copy_array.Rh = self.Rh
        copy_array.Rv = self.Rv
        copy_array.R = np.copy(self.R)
        copy_array.temperature = self.temperature
        copy_array.rightWorkLen = self.rightWorkLen
        copy_array.downWorkLen = self.downWorkLen
        copy_array.no_tunneling_next_time = self.no_tunneling_next_time
        copy_array.voltageToChargeMatrix = self.voltageToChargeMatrix
        copy_array.chargeToVoltageMatrix = self.chargeToVoltageMatrix
        copy_array.equilibriumVoltageMatrix = self.equilibriumVoltageMatrix
        copy_array.invCmat = self.invCmat
        copy_array.chargingEigenVals = self.chargingEigenVals
        copy_array.chargingEigenVecs = self.chargingEigenVecs
        copy_array.invChargingEigenVecs = self.invChargingEigenVecs
        copy_array.deltaV = self.deltaV
        copy_array.deltaVDiag = self.deltaVDiag
        copy_array.deltaVequilibrium = self.deltaVequilibrium
        copy_array.work = np.copy(self.work)
        copy_array.timeStep = self.timeStep
        copy_array.default_dt = self.default_dt
        copy_array.tauLeaping = self.tauLeaping
        copy_array.use_modifyR = self.use_modifyR
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

    def changeVext(self, newVL, newVR, newVU, newVD):
        Q = self.getGroundCharge()
        self.VL = newVL
        self.VR = newVR
        self.VU = newVU
        self.VD = newVD
        self.setInitialVoltages(Q)
        return True

    def getOccupation(self):
        return np.copy(self.n).reshape((self.rows, self.columns))

    def getGroundCharge(self):
        Q = self.voltageToChargeMatrix.dot(self.voltageDiag) - self.getNprime()
        return Q.reshape((self.rows, self.columns))

    def setOccupation(self, n):
        self.n = flattenToColumn(n)
        self.setEquilibriumVoltages()

    def setGroundCharge(self, Q):
        self.voltageDiag = self.chargeToVoltageMatrix.dot(flattenToColumn(Q) + self.getNprime())

    def getTimeStep(self):
        return self.timeStep

    def getCurrentMap(self):
        Ih = self.rates[:self.rightWorkLen] - \
                       self.rates[self.rightWorkLen:2 * self.rightWorkLen]
        Iv = self.rates[2 * self.rightWorkLen:(2 * self.rightWorkLen + self.downWorkLen)] - \
                      self.rates[(2 * self.rightWorkLen + self.downWorkLen):]
        return Ih.reshape((self.rows, self.columns + 1)), Iv.reshape((self.rows + 1, self.columns))

    def getCapacitanceMatrix(self):
        """
        Creates the inverse capacitance matrix
        """
        diagonal = self.Ch[:,:-1] + self.Ch[:,1:] + self.Cv[:-1,:] + self.Cv[1:,:]
        second_diagonal = np.copy(self.Ch[:,1:])
        second_diagonal[:,-1] = 0
        second_diagonal = second_diagonal.flatten()
        second_diagonal = second_diagonal[:-1]
        n_diagonal = np.copy(self.Cv[1:-1,:])
        C_mat = np.diagflat(diagonal) - np.diagflat(second_diagonal,k=1) - np.diagflat(second_diagonal,k=-1)\
                -np.diagflat(n_diagonal, k=self.columns) - np.diagflat(n_diagonal, k=-self.columns)
        return C_mat

    def createConversionMatrices(self):
        cMat = self.getCapacitanceMatrix()
        self.invCmat = np.linalg.inv(cMat)
        X = np.arange(self.columns)
        Y = np.arange(self.rows)
        X, Y = np.meshgrid(X, Y)
        v = self.invCmat[:,5050]
        Z = v.reshape(self.rows, self.columns)
        x = np.arange(self.columns//2)
        y = Z[50:,50]
        # Z2 = np.array([[1/np.sqrt(((i-50.1)**2 + (j-50.1)**2)) for j in range(self.columns)] for i in range(self.rows)])
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(X, Y, Z)
        # ax.plot_surface(X, Y, Z2)
        fig = plt.figure()
        plt.semilogx(x,y)
        plt.show()
        equilibriumChrageMatrix, chargingRateMatrix = self.getChargingRatesmatrix(self.invCmat)
        self.chargingEigenVals, self.chargingEigenVecs, self.invChargingEigenVecs = self.setDiagonalize(self.invCmat.dot(chargingRateMatrix.dot(cMat)))
        self.voltageToChargeMatrix = cMat.dot(self.chargingEigenVecs)
        self.chargeToVoltageMatrix = np.linalg.inv(self.voltageToChargeMatrix)
        self.equilibriumVoltageMatrix = self.invChargingEigenVecs.dot(self.invCmat.dot(np.linalg.inv(equilibriumChrageMatrix)))
        self.timeStep, self.default_dt = self.getTimeConstants(self.chargingEigenVals)

    def getNprime(self):
        left_part_n_prime = np.copy(self.Ch[:, :-1])
        left_part_n_prime[:, 1:] = 0
        right_part_n_prime = np.copy(self.Ch[:, 1:])
        right_part_n_prime[:, :-1] = 0
        up_part_n_prime = np.copy(self.Cv[:-1, :])
        up_part_n_prime[1:, :] = 0
        down_part_n_prime = np.copy(self.Cv[1:, :])
        down_part_n_prime[:-1, :] = 0
        return self.n + flattenToColumn(left_part_n_prime*self.VL +
               right_part_n_prime*self.VR +
               up_part_n_prime*self.VU +
               down_part_n_prime*self.VD)

    def getChargingRatesmatrix(self, invCmat):
        """
        For developement of Q by dQ/dt=JQ +b
        """
        res = invCmat + np.diagflat(1/self.CG)
        return res, res / np.repeat(flattenToColumn(self.RG),res.shape[1],axis=1)

    def setDiagonalize(self, mat):
        eigenValues, eigenVectors = np.linalg.eig(mat)
        eigenValues = flattenToColumn(eigenValues)
        eigenVectorsInv = np.linalg.inv(eigenVectors)
        return eigenValues, eigenVectors, eigenVectorsInv

    def getTimeConstants(self, chargingRateEigenValues):
        timeStep = 1/np.min(chargingRateEigenValues)
        default_dt = 0.1/np.max(chargingRateEigenValues)
        return timeStep, default_dt

    def setInitialVoltages(self, Q0):
        self.voltages = self.invCmat.dot(self.getNprime() + flattenToColumn(Q0))
        self.voltagesDiag = self.invChargingEigenVecs.dot(self.voltages)

    def developeV(self, dt):
        self.voltageDiag = (self.voltagesDiag - self.equilibriumVoltages) * np.exp(self.chargingEigenVals) +\
                           self.equilibriumVoltages
        return True

    def getVoltagesDiag(self):
        return self.voltageDiag

    def loadState(self, n, voltages):
        self.n = np.copy(n)
        self.voltages = np.copy(voltages)
        self.voltageDiag = self.invChargingEigenVecs.dot(voltages)
        self.setEquilibriumVoltages()


    def setVoltages(self):
        self.voltages = self.chargingEigenVecs.dot(self.voltagesDiag)

    def getVoltages(self):
        return self.voltages.reshape((self.rows,self.columns))

    def getDistFromEquilibrium(self):
        return self.voltageDiag - self.equilibriumVoltages

    def convertVdiagToV(self, vDiag):
        return self.chargingEigenVecs.dot(vDiag)

    def setEquilibriumVoltages(self):
        self.equilibriumVoltages = self.equilibriumVoltageMatrix.dot(self.VG + self.getNprime()/self.CG)

    def setChangeMatrices(self):
        dVright, dVdown = self.changeHelper(self.invCmat)
        self.deltaV = np.hstack((dVright,-dVright,dVdown,-dVdown))
        dVDiagright, dVDiagdown = self.changeHelper(self.chargeToVoltageMatrix)
        self.deltaVDiag = np.hstack((dVDiagright,-dVDiagright,dVDiagdown,-dVDiagdown))
        dVequilibriumRight, dVequilibriumDown = self.changeHelper(self.equilibriumVoltageMatrix.dot(np.diagflat(1/self.CG)))
        self.deltaVequilibrium = np.hstack((dVequilibriumRight, -dVequilibriumRight, dVequilibriumDown, -dVequilibriumDown))

    def changeHelper(self, mat):
        padded = np.pad(mat.reshape((self.rows,self.columns,self.rows,self.columns)), ((0,0),(0,0),(1,1),(1,1)), mode="constant")
        right = padded[:,:,1:-1,1:] - padded[:,:,1:-1,:-1]
        down = padded[:,:,1:,1:-1] - padded[:,:,:-1,1:-1]
        return right.reshape((self.rows*self.columns,self.rightWorkLen)), down.reshape((self.rows*self.columns,self.downWorkLen))

    def getWork(self, charge=1):
        padded = np.pad(self.voltages.reshape((self.rows, self.columns)), ((1,1),(1,1)), mode='constant',
                        constant_values=((self.VL,self.VR),(self.VU,self.VD)))
        Vright = (padded[1:-1,1:] - padded[1:-1,:-1]).flatten()
        Vdown = (padded[1:,1:-1] - padded[:-1,1:-1]).flatten()
        dVright = (np.pad(np.diagonal(self.deltaV[:,:self.rows*self.columns]).reshape((self.rows, self.columns)),((0,0), (0,1))) - \
                  np.pad(np.diagonal(self.deltaV[:,1:self.rows*self.columns+1]).reshape((self.rows, self.columns)),((0,0),(1,0)))).flatten()
        dVdown = (np.pad(np.diagonal(self.deltaV[:,self.rightWorkLen:self.rightWorkLen+self.rows*self.columns]).reshape((self.rows, self.columns)),((0,1),(0,0))) - \
                 np.pad(np.diagonal(self.deltaV[:,self.rightWorkLen+1:self.rightWorkLen+1+self.rows*self.columns]).reshape((self.rows, self.columns)),((1,0),(0,0)))).flatten()
        self.work[:self.rightWorkLen] = charge*Vright + 0.5*(charge**2)*dVright
        self.work[self.rightWorkLen:2*self.rightWorkLen] = -charge*Vright + 0.5*(charge**2)*dVright
        self.work[2*self.rightWorkLen:2*self.rightWorkLen+self.downWorkLen] = charge*Vdown + 0.5*(charge**2)*dVdown
        self.work[2 * self.rightWorkLen + self.downWorkLen:] = -charge*Vdown + 0.5*(charge**2)*dVdown
        return self.work

    def updateV(self, actionId, charge):
        self.voltages = self.voltages + flattenToColumn(charge*self.deltaV[:,actionId])
        self.voltagesDiag = self.voltagesDiag + flattenToColumn(charge*self.deltaVDiag[:,actionId])
        self.equilibriumVoltages = self.equilibriumVoltages + flattenToColumn(charge * self.deltaVequilibrium[:, actionId])

    def updateVMultipleActions(self, actionVec, charge):
        deltaVDiag = np.zeros(self.voltagesDiag.shape)
        deltaVequilibrium = np.zeros(self.equilibriumVoltages.shape)
        for ind, times in enumerate(actionVec):
            deltaVDiag += times * self.deltaVDiag[:, ind]
            deltaVequilibrium += times * self.deltaVequilibrium[:, ind]
        self.voltagesDiag += charge * deltaVDiag
        self.equilibriumVoltages += charge * deltaVequilibrium

    def getRates(self):
        """
        Returns the tunnelling rate between neighboring dots
        :param T: temperature (in energy units = e^2[V])
        """
        work = self.getWork()
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
        rightCurrent = np.sum(self.rates[:self.rightWorkLen]) -\
                      np.sum(self.rates[self.rightWorkLen:2*self.rightWorkLen])
        downCurrent = np.sum(self.rates[2*self.rightWorkLen:(2*self.rightWorkLen + self.downWorkLen)]) -\
                      np.sum(self.rates[(2*self.rightWorkLen + self.downWorkLen):])
        return rightCurrent/(self.columns+1), downCurrent/(self.rows+1)

    def modifyR(self):
        nExponent = np.exp(-INV_DOS*self.n).reshape((self.rows, self.columns))
        self._rightnExponent[:,1:] = nExponent
        self._leftnExponent[:,:-1] = nExponent
        horzSize = self.Rh.size
        vertSize = self.Rv.size
        self.R[:horzSize] = (self.Rh*self._rightnExponent).flatten()
        self.R[horzSize:2*horzSize] = (self.Rh * self._leftnExponent).flatten()
        if vertSize > 0:
            self.R[2*horzSize:2*horzSize+vertSize] = (self.Rv*nExponent[:-1,:]).flatten()
            self.R[2 * horzSize + vertSize:] = (self.Rv * nExponent[1:, :]).flatten()


    def getVoltagesFromGround(self):
        return self.VG.reshape(self.rows, self.columns) - self.getGroundCharge()/self.CG.reshape(self.rows, self.columns)

    def get_dist_from_steady(self):
        return np.abs(self.qDiag - self.qDiagEq)

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
        out_of_interval = dt > self.default_dt
        while out_of_interval:
            intervals += 1
            self.developeV(self.default_dt)
            self.setVoltages()
            rates = self.getRates()
            np.cumsum(rates, out=self.cumRates)
            new_sum_rates = self.cumRates[-1]
            if new_sum_rates == 0:
                self.no_tunneling_next_time = True
                return self.default_dt, intervals
            dt = (dt - self.default_dt)*(new_sum_rates/sum_rates)
            sum_rates = new_sum_rates
            out_of_interval = dt > self.default_dt
            if np.max(self.get_dist_from_steady()) < ALLOWED_ERR:
                break
        return dt, intervals

    def getLeapingTimeInterval(self):
        rates = self.getRates()
        if (rates <= EPS).all():
            self.no_tunneling_next_time = True
            return self.default_dt
        horzSize = self.rows * (self.columns + 1)
        vertSize = (self.rows + 1) * self.columns
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
        self.developeV(dt)
        self.setVoltages()
        if self.no_tunneling_next_time:
            self.no_tunneling_next_time = False
            return True
        if (self.rates == 0).all():
            return True
        actionInd = np.searchsorted(self.cumRates, randomNumber*self.cumRates[-1])
        self.executeAction(actionInd)
        return True

    def nextLeapingSteps(self, dt):
        self.developeV(dt)
        self.setVoltages()
        if self.no_tunneling_next_time:
            self.no_tunneling_next_time = False
            return True
        actionVec = np.random.poisson(lam=self.rates*dt, size=self.rates.size)
        self.executeMultipleActions(actionVec)

    def tunnel(self, fromDot, toDot, charge=1):
        if (0 <= fromDot[1] < self.columns) and (0 <= fromDot[0] < self.rows):
            self.n[fromDot[0]*self.columns + fromDot[1],0] -= charge
        if (0 <= toDot[1] < self.columns) and (0 <= toDot[0] < self.rows):
            self.n[toDot[0]*self.columns + toDot[1],0] += charge
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
        vertSize = (self.rows + 1) * self.columns
        self.right[:,:] = actionVec[:horzSize].reshape(self.right.shape)*charge
        self.left[:,:] = actionVec[horzSize:2 * horzSize].reshape(self.left.shape)*charge
        self.down[:,:] = actionVec[2 * horzSize:2 * horzSize + vertSize].reshape(self.down.shape)*charge
        self.up[:,:] = actionVec[2 * horzSize + vertSize:].reshape(self.up.shape)*charge
        self.totalAction[:,:] = 0
        self.totalAction += self.left[:,1:] - self.left[:,:-1] + self.right[:,:-1] - self.right[:,1:]
        self.totalAction[1:,:] += self.down - self.up
        self.totalAction[:-1,:] += self.up - self.down
        self.n += flattenToColumn(self.totalAction)
        self.updateVMultipleActions(actionVec, charge)

    def executeAction(self, ind, charge=1):
        horzSize = self.rows*(self.columns+1)
        vertSize = (self.rows + 1)*self.columns
        if ind < horzSize: # tunnel right:
            fromDot = (ind//(self.columns+1), (ind%(self.columns+1))-1)
            toDot = (fromDot[0], fromDot[1] + 1)
        elif ind < horzSize*2: # tunnel left
            ind -= horzSize
            fromDot = (ind//(self.columns+1), ind%(self.columns+1))
            toDot = (fromDot[0], fromDot[1] - 1)
        elif ind < horzSize*2 + vertSize: # tunnel down
            ind -= horzSize*2
            fromDot = (ind//self.columns-1, ind%self.columns)
            toDot = (fromDot[0] + 1, fromDot[1])
        else: # tunnel up
            ind -= (horzSize*2 + vertSize)
            fromDot = ((ind//self.columns), ind%self.columns)
            toDot = (fromDot[0] - 1, fromDot[1])
        self.tunnel(fromDot, toDot, charge=charge)
        self.updateV(ind, charge)
        return fromDot, toDot

class JJArray(DotArray):
    def __init__(self, rows, columns, VL, VR, VU, VD,  VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv,
                 temperature, scGap, fastRelaxation=False, tauLeaping=False, modifyR=False):
        DotArray.__init__(self, rows, columns, VL, VR, VU, VD, VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv,
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

    def getWork(self, charge=1):
        qp_work = super().getWork(charge=1)
        cp_work = super().getWork(charge=2)
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
        vertSize = (self.rows + 1)*self.columns
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
    def __init__(self, index, VL0, VR0, VU0, VD0, Q0, n0, dotArray, constQ=False):
        self.dotArray = copy(dotArray)
        self.randomGenerator = np.random.RandomState()
        self.dotArray.setGroundCharge(Q0)
        self.dotArray.setOccupation(n0)
        self.VL = VL0
        self.VR = VR0
        self.VD = VD0
        self.VU = VU0
        self.index = index
        self.tauLeaping = dotArray.tauLeaping
        self.constQ = constQ
        self.minSteps = MIN_STEPS*self.dotArray.columns
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
        if self.constQ:
            self.find_next_QG_using_gradient_descent()
            self.dotArray.setOccupation(self.n)
            self.dotArray.setGroundCharge(self.Q)
            self.dotArray.getWork()
        final_t = self.dotArray.timeStep
        curr_t = 0
        steps = 0
        I_avg = 0
        I_var = 0
        vert_I_avg = 0
        vert_I_var = 0
        if fullOutput:
            n_avg = np.zeros((self.dotArray.getRows(), self.dotArray.getColumns()))
            n_var = np.zeros(n_avg.shape)
            V_avg = np.zeros((self.dotArray.getRows(), self.dotArray.getColumns()))
            V_var = np.zeros(V_avg.shape)
        if currentMap:
            Ih_avg = np.zeros((self.dotArray.getRows(), self.dotArray.getColumns() + 1))
            Iv_avg = np.zeros((self.dotArray.getRows() + 1, self.dotArray.getColumns()))
        curr_n = self.dotArray.getOccupation()
        curr_V = self.dotArray.getVoltages()
        # err = 2*ALLOWED_ERR
        # errors = []
        while curr_t < final_t or steps < MIN_STEPS:
            if self.tauLeaping:
                dt = self.executeLeapingStep(printState=print_stats)
            else:
                dt = self.executeStep(printState=print_stats)
            steps += 1
            curr_I, curr_vert_I = self.dotArray.getCurrentFromRates()
            I_avg, I_var = self.update_statistics(curr_I, I_avg, I_var,
                                                  curr_t, dt)
            vert_I_avg, vert_I_var = self.update_statistics(curr_vert_I, vert_I_avg, vert_I_var,
                                                  curr_t, dt)
            if fullOutput:
                n_avg, n_var = self.update_statistics(curr_n, n_avg, n_var, curr_t, dt)
                V_avg, V_var = self.update_statistics(curr_V, V_avg, V_var, curr_t, dt)
                curr_n = self.dotArray.getOccupation()
                curr_Q = self.dotArray.getGroundCharge()
            if currentMap:
                Ih, Iv = self.dotArray.getCurrentMap()
                Ih_avg += Ih*dt
                Iv_avg += Iv*dt
            curr_t += dt
        err = self.get_err(I_var, steps, curr_t)
        vert_err = self.get_err(vert_I_var, steps, curr_t)
        res = (I_avg, err, vert_I_avg, vert_err)
        if fullOutput:
            V_avg = self.dotArray.convertVdiagToV(V_avg)
            V_err = self.dotArray.convertVdiagToV(self.get_err(V_var, steps, curr_t))
            n_err = self.get_err(n_var, steps, curr_t)
            res = res + (n_avg, V_avg, n_err, V_err)
        if currentMap:
            res = res + (Ih/curr_t, Iv/curr_t)
        return res

    def getToSteadyState(self):
        curr_t = 0
        steps = 0
        err = ALLOWED_ERR*2
        not_decreasing = 0
        curr_dist = self.dotArray.getDistFromEquilibrium()
        dist_avg = 0
        dist_var = 0
        # plot = True
        # if plot:
        #     distances = []
        #     ts = []
        #     t=0
        while err > ALLOWED_ERR and not_decreasing < STEADY_STATE_REP:
            if self.tauLeaping:
                dt = self.executeLeapingStep()
            else:
                dt = self.executeStep()
            steps += 1
            dist_avg, dist_var = self.update_statistics(curr_dist, dist_avg, dist_var,
                                                  curr_t, dt)
            curr_dist = self.dotArray.getDistFromEquilibrium()
            curr_t += dt
            # if plot:
            #     t+=dt
            #     distances.append(curr_dist)
            #     ts.append(t)
            if steps % MIN_STEPS == 0:
                new_err = np.max(dist_avg)
                if err < new_err:
                    not_decreasing += 1
                err = new_err
                n_avg = np.zeros(
                    (self.dotArray.getRows(), self.dotArray.getColumns()))
                dist_avg = 0
                dist_var = 0
                curr_t = 0
                # if plot and steps % (100 * MIN_STEPS):
                #     distances = np.array(diatances)
                #     for i in range(distances.shape[2]):
                #         plt.plot(ts, distances[:,:,i],'.')
                #     print(err)
                #     plt.show()
                #     distances = []
                #     ts = []

        return True

    # def checkSteadyState(self, rep):
    #     I = []
    #     ns = []
    #     Qs = []
    #     t = self.dotArray.getTimeStep()/10
    #     for r in range(rep):
    #         # running once to get to steady state
    #         self.calcCurrent(t)
    #         # now we are in steady state calculate current
    #         current, currentErr, n, Q, nErr, QErr =\
    #             self.calcCurrent(t, fullOutput=True)
    #         I.append(current)
    #         ns.append(n)
    #         Qs.append(Q)
    #     ns = np.array(ns)
    #     Qs = np.array(Qs)
    #     x = t*np.arange(rep)
    #     print("Mean: " + str(np.average(I))+ " Var right: " + str(np.var(I)))
    #     print("Tail mean : " + str(np.average(I[len(I)//4:])) + " Tail var right: " + str(np.var(I[len(I)//4:])))
    #     plt.figure()
    #     plt.plot(x,I)
    #     plt.figure()
    #     plt.plot(x, ns.reshape((ns.shape[0], ns.shape[1]*ns.shape[2])))
    #     plt.figure()
    #     plt.plot(x, Qs.reshape((Qs.shape[0], Qs.shape[1]*Qs.shape[2])))
    #     plt.show()

    def saveState(self, I, IErr, vertI, vertIErr, n=None, V=None, nErr=None, VErr=None, Ih=None, Iv=None,
                  fullOutput=False, currentMap=False, basePath=''):
        baseName = basePath + "_temp_" + str(self.index)
        if fullOutput:
            np.save(baseName + "_ns", np.array(n))
            np.save(baseName + "_Vs", np.array(V))
            np.save(baseName + "_nsErr", np.array(nErr))
            np.save(baseName + "_VsErr", np.array(VErr))
        np.save(baseName + "_n", self.dotArray.getOccupation())
        np.save(baseName + "_V", self.dotArray.getVoltages())
        if currentMap:
            np.save(baseName + "_Ih", np.array(Ih))
            np.save(baseName + "_Iv", np.array(Iv))
        np.save(baseName + "_I", np.array(I))
        np.save(baseName + "_IErr", np.array(IErr))
        np.save(baseName + "_vertI", np.array(vertI))
        np.save(baseName + "_vertIErr", np.array(vertIErr))

    def loadState(self,  fullOutput=False, currentMap=False, basePath=''):
        baseName = basePath + "_temp_" + str(self.index)
        I = np.load(baseName + "_I.npy")
        loadLen = len(I)
        IErr = np.load(baseName + "_IErr.npy")
        vertI = np.load(baseName + "_vertI.npy")
        vertIErr = np.load(baseName + "_vertIErr.npy")
        if len(IErr) > loadLen:
            IErr = IErr[:loadLen]
        if len(vertI) > loadLen:
            vertI = vertI[:loadLen]
        if len(vertIErr) > loadLen:
            vertIErr = vertIErr[:loadLen]
        n = np.load(baseName + "_n.npy")
        V = np.load(baseName + "_V.npy")
        res = (I,IErr,vertI,vertIErr,n,V)
        if fullOutput:
            ns = np.load(baseName + "_ns.npy")
            if len(ns) > loadLen:
                ns = ns[loadLen, :, :]
            Vs = np.load(baseName + "_Vs.npy")
            if len(Vs) > loadLen:
                Vs = Vs[loadLen, :, :]
            nsErr = np.load(baseName + "_nsErr.npy")
            if len(nsErr) > loadLen:
                nsErr = nsErr[loadLen, :, :]
            VsErr = np.load(baseName + "_VsErr.npy")
            if len(VsErr) > loadLen:
                VsErr = VsErr[:, loadLen, :, :]
            res = res + (ns, Vs,nsErr, VsErr)
        if currentMap:
            Ih = np.load(baseName + "_Ih.npy")
            Iv = np.load(baseName + "_Iv.npy")
            if len(Ih) > loadLen:
                Ih = Ih[:, loadLen, :, :]
            if len(Iv) > loadLen:
                Iv = Iv[:, loadLen, :, :]
            res = res + (Ih,Iv,)
        return res

    def calcIV(self, Vmax, Vstep, vSym, fullOutput=False, print_stats=False,
               currentMap=False, basePath="", resume=False):
        I = []
        IErr = []
        vertI = []
        vertIErr = []
        ns = []
        voltages = []
        nsErr = []
        voltagesErr = []
        Ih = []
        Iv = []
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
            vertI = list(resumeParams[2])
            vertIErr = list(resumeParams[3])
            Vind = len(I)
            VL_vec = VL_vec[Vind:]
            VR_vec = VR_vec[Vind:]
            occupation = resumeParams[4]
            vs = resumeParams[5]
            self.dotArray.loadState(occupation, vs)
            if fullOutput:
                ns = list(resumeParams[6])
                voltages = list(resumeParams[7])
                nsErr = list(resumeParams[8])
                voltagesErr = list(resumeParams[9])
            if currentMap:
                Ih = list(resumeParams[-2])
                Iv = list(resumeParams[-1])
        for VL,VR in zip(VL_vec, VR_vec):
            self.dotArray.changeVext(VL, VR, self.VU, self.VD)
            # running once to get to steady state
            if not self.constQ:
                self.getToSteadyState()
            # now we are in steady state calculate current
            stepRes = self.calcCurrent(print_stats=print_stats, fullOutput=fullOutput, currentMap=currentMap)
            current = stepRes[0]
            currentErr = stepRes[1]
            vert_current = stepRes[2]
            vert_currentErr = stepRes[3]
            if fullOutput:
                ns.append(stepRes[4])
                voltages.append(stepRes[5])
                nsErr.append(stepRes[6])
                voltagesErr.append(stepRes[7])
            if currentMap:
                Ih.append(stepRes[-2])
                Iv.append(stepRes[-1])
            I.append(current)
            IErr.append(currentErr)
            vertI.append(vert_current)
            vertIErr.append(vert_currentErr)
            if self.index == 0:
                print(VL - VR, end=',', flush=True)
            self.saveState(I, IErr, vertI, vertIErr, ns, voltages, nsErr, voltagesErr, Ih=Ih, Iv=Iv, fullOutput=fullOutput,
                           currentMap=currentMap, basePath=basePath)
        res = (np.array(I), np.array(IErr), np.array(vertI), np.array(vertIErr), VL_res - VR_res)
        if fullOutput:
            res = res + (np.array(ns), np.array(voltages), np.array(nsErr), np.array(voltagesErr))
        if currentMap:
            res = res + (np.array(Ih),np.array(Iv))
        return res

    def printState(self):
        self.dotArray.printState()



class GraphSimulator:
    """
    Calculating steady state current for an array of quantum dots by using linear programming
    """
    def __init__(self, index, VL0, VR0, VU0, VD0, Q0, n0, dotArray):
        self.dotArray = copy(dotArray)
        self.n0 = n0
        self.QG = Q0
        self.VL = VL0
        self.VR = VR0
        self.VU = VU0
        self.VD= VD0
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
            self.dotArray.changeVext(VL, VR, self.VU, self.VD)
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

def runSingleSimulation(index, VL0, VR0, VU0, VD0, vSym, Q0, n0,Vmax, Vstep, dotArray,
                        fullOutput=False, printState=False, useGraph=False, currentMap=False,
                        basePath="", resume=False, constQ=False):
    if useGraph:
        simulator = GraphSimulator(index, VL0, VR0, VU0, VD0, Q0, n0, dotArray)
    else:
        simulator = Simulator(index, VL0, VR0, VU0, VD0, Q0, n0, dotArray, constQ)
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

def runFullSimulation(VL0, VR0, VU0, VD0, vSym, VG0, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, rows, columns,
                      Vmax, Vstep, temperature=0, repeats=1, savePath=".", fileName="", fullOutput=False,
                      printState=False, checkSteadyState=False, useGraph=False, fastRelaxation=False,
                      currentMap=False, dbg=False, plotCurrentMaps=False, plotBinaryCurrentMaps=False, resume=False, superconducting=False,
                      gap=0, leaping=False, modifyR=False, plotVoltages=False,
                      constQ=False):
    basePath = os.path.join(savePath, fileName)
    if useGraph:
        fastRelaxation = True
        repeats = 1
    if plotCurrentMaps or plotBinaryCurrentMaps:
        print("Plotting Current Maps")
        Ih = np.load(basePath + "_Ih.npy")
        Iv = np.load(basePath + "_Iv.npy")
        V = np.load(basePath + "_V.npy")
        n=None
        if fullOutput:
            n = np.load(basePath + "_n.npy")
        saveCurrentMaps(Ih, Iv, V, basePath + "_Imap",full=fullOutput,
                        n=n, binary=plotBinaryCurrentMaps)
        exit(0)
    # if not resume:
    #     print("Saving array parameters")
    #     saveRandomParams(np.array(VG0),
    #                      np.array(Q0), np.array(n0), np.array(CG),
    #                      np.array(RG), np.array(Ch), np.array(Cv),
    #                      np.array(Rh), np.array(Rv), basePath)
    # else:
    #     print("Loading array parameters")
    #     VG0, Q0, n0, CG, RG, Ch, Cv, Rh, Rv = loadRandomParams(basePath)
    if superconducting:
        prototypeArray = JJArray(rows, columns, VL0, VR0, VU0, VD0, np.array(VG0), np.array(Q0), np.array(n0), np.array(CG),
                                 np.array(RG), np.array(Ch), np.array(Cv), np.array(Rh), np.array(Rv),
                                 temperature, gap, fastRelaxation=fastRelaxation, tauLeaping=leaping, modifyR=modifyR)
        print("Superconducting prototype array was created")
    else:
        prototypeArray = DotArray(rows, columns, VL0, VR0, VU0, VD0, np.array(VG0), np.array(Q0), np.array(n0), np.array(CG),
                                  np.array(RG), np.array(Ch), np.array(Cv), np.array(Rh), np.array(Rv),
                                  temperature, fastRelaxation=fastRelaxation, tauLeaping=leaping, modifyR=modifyR,
                                  constQ=constQ)
        print("Normal prototype array was created")
    if checkSteadyState:
        # print("Running steady state convergance check")
        # simulator = Simulator(0, VL0, VR0, VU0, VD0, Q0, n0, prototypeArray)
        # simulator.checkSteadyState(repeats)
        # simulator.dotArray.changeVext(VL0+Vstep, VR0, VU0, VD0)
        # simulator.checkSteadyState(repeats)
        exit(0)

    # if plotVoltages:
        # print("Plotting voltages")
        # points_num = 100
        # Q0 = np.array([[1,1,0.1]])
        # dQ = np.array([[0,0,0.01]])
        # VL0 = 1
        # simulator = Simulator(0, VL0, VR0, VU0, VD0, Q0, n0, prototypeArray)
        # simulator.plotAverageVoltages(Q0,dQ,n0,points_num)
        # exit(0)
    Is = []
    IsErr = []
    vertIs = []
    vertIsErr = []
    ns = []
    Vs = []
    nsErr = []
    VsErr = []
    Ih = []
    Iv = []
    if not dbg:
        print("Starting parallel run")
        pool = Pool(processes=repeats)
        results = []
        for repeat in range(repeats):
            res = pool.apply_async(runSingleSimulation,
                                    (repeat, VL0, VR0, VU0, VD0, vSym, Q0, n0, Vmax, Vstep, prototypeArray, fullOutput,
                                     printState, useGraph, currentMap,basePath, resume, constQ))
            results.append(res)
        for res in results:
            result = res.get()
            I = result[0]
            IErr = result[1]
            vertI = result[2]
            vertIErr = result[3]
            if fullOutput:
                n = result[5]
                voltages = result[6]
                nErr = result[7]
                VErr = result[8]
                ns.append(n)
                Vs.append(voltages)
                nsErr.append(nErr)
                VsErr.append(VErr)
            if currentMap:
                Ih.append(result[-3])
                Iv.append(result[-2])
            Is.append(I)
            IsErr.append(IErr)
            vertIs.append(vertI)
            vertIsErr.append(vertIErr)
        V = result[4]
        params = result[-1]

    else: #dbg
        print("Starting serial run")
        for repeat in range(repeats):
            result = runSingleSimulation(repeat, VL0, VR0, VU0, VD0, vSym, Q0, n0, Vmax, Vstep, prototypeArray, fullOutput,
                                         printState, useGraph, currentMap,basePath, resume, constQ)
            I = result[0]
            IErr = result[1]
            vertI = result[2]
            vertIErr = result[3]
            if fullOutput:
                n = result[5]
                voltages = result[6]
                nErr = result[7]
                VErr = result[8]
                ns.append(n)
                Vs.append(voltages)
                nsErr.append(nErr)
                VsErr.append(VErr)
            if currentMap:
                Ih.append(result[-3])
                Iv.append(result[-2])
            Is.append(I)
            IsErr.append(IErr)
            vertIs.append(vertI)
            vertIsErr.append(vertIErr)
        V = result[4]
        params = result[-1]
        # dbg
    print("Saving results")
    avgI = np.mean(np.array(Is), axis=0)
    avgVertI = np.mean(np.array(vertIs), axis=0)
    if repeats < 10:
        avgIErr = np.sqrt(np.sum(np.array(IsErr)**2, axis=0)/len(IsErr))
        avgVertIErr = np.sqrt(np.sum(np.array(vertIsErr)**2, axis=0)/len(vertIsErr))
    else:
        avgIErr = np.std(np.array(Is), axis=0)/np.sqrt(len(Is))
        avgVertIErr = np.std(np.array(vertIs), axis=0) / np.sqrt(len(vertIs))
    if fullOutput:
        avgN = np.mean(np.array(ns), axis=0)
        avgV = np.mean(np.array(Vs), axis=0)
        avgNErr = np.sqrt(np.sum(np.array(nsErr) ** 2, axis=0)) / len(nsErr)
        avgVErr = np.sqrt(np.sum(np.array(VsErr) ** 2, axis=0)) / len(QsErr)
        np.save(basePath + "_n", avgN)
        np.save(basePath + "_V", avgV)
        np.save(basePath + "_nErr", avgNErr)
        np.save(basePath + "_VErr", avgVErr)
        np.save(basePath + "_full_I", np.array(Is))
        np.save(basePath + "_full_IErr", np.array(IsErr))
    fig = plt.figure()
    IplusErr = avgI + avgIErr
    IminusErr = avgI - avgIErr
    plt.plot(V[:V.size//2], avgI[:V.size//2], 'b-',V[:V.size//2], IplusErr[:V.size//2], 'b--',
             V[:V.size//2], IminusErr[:V.size//2], 'b--',
             V[V.size//2:], avgI[V.size//2:], 'r-',V[V.size//2:], IplusErr[V.size//2:], 'r--',
             V[V.size//2:], IminusErr[V.size//2:], 'r--')
    plt.xlabel('Voltage')
    plt.ylabel('Current from left to right')
    plt.savefig(basePath + "_IV.png")
    plt.close(fig)
    fig = plt.figure()
    IplusErr = avgVertI + avgVertIErr
    IminusErr = avgVertI - avgVertIErr
    plt.plot(V[:V.size // 2], avgVertI[:V.size // 2], 'b-', V[:V.size // 2],
             IplusErr[:V.size // 2], 'b--',
             V[:V.size // 2], IminusErr[:V.size // 2], 'b--',
             V[V.size // 2:], avgVertI[V.size // 2:], 'r-', V[V.size // 2:],
             IplusErr[V.size // 2:], 'r--',
             V[V.size // 2:], IminusErr[V.size // 2:], 'r--')
    plt.xlabel('Voltage')
    plt.ylabel('Current from up to down')
    plt.savefig(basePath + "_vert_IV.png")
    plt.close(fig)
    np.save(basePath + "_I", avgI)
    np.save(basePath + "_IErr", avgIErr)
    np.save(basePath + "_vertI", avgVertI)
    np.save(basePath + "_vertIErr", avgVertIErr)
    np.save(basePath + "_V", V)

    if currentMap:
        avgIh = np.mean(np.array(Ih),axis=0) if repeats > 1 else np.array(Ih[0])
        np.save(basePath + "_Ih", avgIh)
        avgIv = np.mean(np.array(Iv), axis=0) if repeats > 1 else np.array(Iv[0])
        np.save(basePath + "_Iv", avgIv)
    for index in range(repeats):
        removeState(index, fullOutput=fullOutput, basePath=basePath, currentMap=currentMap, graph=use_graph)
    removeRandomParams(basePath)
    return params

def removeState(index, fullOutput=False, basePath='', currentMap=False, graph=False):
    baseName = basePath + "_temp_" + str(index)
    os.remove(baseName + "_I.npy")
    if not graph:
        os.remove(baseName + "_IErr.npy")
        os.remove(baseName + "_vertI.npy")
        os.remove(baseName + "_vertIErr.npy")
    os.remove(baseName + "_n.npy")
    os.remove(baseName + "_V.npy")
    if fullOutput:
        os.remove(baseName + "_ns.npy")
        os.remove(baseName + "_Vs.npy")
        if not graph:
            os.remove(baseName + "_nsErr.npy")
            os.remove(baseName + "_VsErr.npy")
    if currentMap:
        os.remove(baseName + "_Ih.npy")
        os.remove(baseName + "_Iv.npy")
    return True

def saveCurrentMaps(Ih, Iv, V, path, full=False, n=None, binary=False):
    if binary:
        Ih[Ih>0] = 1
        Ih[Ih<0] = -1
        Iv[Iv > 0] = 1
        Iv[Iv < 0] = -1
        path += "_binary"
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=24, bitrate=1800)
    fig, ax = plt.subplots()
    Imax = max(np.max(Ih), np.max(Iv))
    M,N = Ih[0].shape
    N -= 1
    im = ax.imshow(np.zeros((M*3 + 2, N * 3+2)), vmin=-Imax,
                    vmax=Imax, animated=True, cmap='PuOr', aspect='equal')
    text = ax.text(1, 1, 'Vext = 0')
    cb1 = plt.colorbar(im, shrink=0.25)
    cb1.set_label('Current')
    if full:
        nmax = np.max(n)
        nmin = np.min(n)
        im2 = ax.imshow(np.zeros((M*3 + 2, N * 3+2)),
                        vmin=nmin, vmax=nmax, animated=True, cmap='RdBu',
                        aspect='equal')
        cb2 = plt.colorbar(im2, shrink=0.25)
        cb2.set_label('Occupation')
        frames = [(Ih[i], Iv[i], V[i], n[i]) for i in range(len(V))]
        im_ani = animation.FuncAnimation(fig,
                                         plotCurrentMaps(im, text, M, N, full=True, im2=im2),
                                         frames=frames, interval=100,
                                         repeat_delay=1000,
                                         blit=True)
    else:
        frames = [(Ih[i], Iv[i], V[i]) for i in range(len(V))]
        im_ani = animation.FuncAnimation(fig, plotCurrentMaps(im, text, M, N),
                                        frames=frames, interval=100,
                                         repeat_delay=1000,
                                        blit=True)
    im_ani.save(path + '.mp4', writer=writer)
    plt.close(fig)

def plotCurrentMaps(im, text, M, N, full=False, im2=None):
    '''
    updating the plot to current currents map
    :return: image for animation
    '''
    J = np.zeros((M*3+2,N*3+2))
    horzRows = np.arange(2,M*3,3)
    horzCols = np.repeat(np.arange(0,3*N+1,3),2)
    horzCols[1::2] += 1
    vertRows = np.repeat(np.arange(0, M*3+1, 3),2)
    vertRows[1::2] += 1
    vertCols = np.arange(2,3*N,3)
    Jmask = np.ones(J.shape)
    Jmask[np.ix_(horzRows, horzCols)] = 0
    Jmask[np.ix_(vertRows, vertCols)] = 0
    if full:
        dot_rows = horzRows
        dot_cols = vertCols
        dots_im = J.copy()
        dots_im_mask = np.ones(J.shape)
        dots_im_mask[np.ix_(dot_rows, dot_cols)] = 0
    def updateCurrent(result):
        if full:
            Ih, Iv, Vext, n = result
        else:
            Ih, Iv, Vext = result
        if Ih is None:
            return im
        J[np.ix_(horzRows, horzCols)] = np.repeat(Ih,2,axis=1)
        J[np.ix_(vertRows, vertCols)] = np.repeat(Iv,2,axis=0)
        J_masked = np.ma.masked_array(J,Jmask)
        im.set_array(J_masked)
        text.set_text('Vext = ' + str(Vext))
        if full:
            dots_im[np.ix_(dot_rows, dot_cols)] = n
            dots_im_masked = np.ma.masked_array(dots_im, dots_im_mask)
            im2.set_array(dots_im_masked)
            return im, text, im2
        else:
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
    parser.add_option("--vu", dest="VU",
                      help="upper electrode voltage (in units of"
                           " planckConstant/electronCharge*timeUnits) [default: %default]",
                      default=0, type=float)
    parser.add_option("--vd", dest="VD",
                      help="lower electrode voltage (in units of"
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
    parser.add_option("--plot-binary-current-map", dest="plot_binary_current_map",
                      help="if true a binary clip of current distribution will"
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
    VU0 = options.VU
    VD0 = options.VD
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
    plot_binary_current_map = options.plot_binary_current_map
    modifyR = options.modifyR
    constQ = options.constQ
    if params_file:
        VG = arrayParams['VG']
        CG = arrayParams['CG']
        RG = arrayParams['RG']
        Ch = arrayParams['Ch']
        Cv = arrayParams['Cv']
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
            Cv = create_random_array(rows + 1, columns, options.C_avg, options.C_std, dist,
                                     True)
        if options.custom_rh.replace('\"', '') and options.custom_rv.replace('\"', ''):
            Rh = literal_eval(options.custom_rh.replace('\"', ''))
            Rv = literal_eval(options.custom_rv.replace('\"', ''))
        else:
            Rh = create_random_array(rows, columns + 1, options.R_avg, options.R_std, dist,
                                     True)
            Rv = create_random_array(rows + 1, columns, options.R_avg, options.R_std, dist,
                                     True)

    # Running Simulation
    if not os.path.exists(savePath):
        os.mkdir(savePath)
    elif not os.path.isdir(savePath):
        print("the given path exists but is a file")
        exit(0)
    if dbg:
        import cProfile, pstats, io
        pr = cProfile.Profile()
        pr.enable()

    array_params = runFullSimulation(VL0, VR0, VU0, VD0, vSym, VG, Q0, n0, CG, RG, Ch, Cv, Rh, Rv, rows,  columns,
                                     Vmax, Vstep, temperature=T, repeats=repeats, savePath=savePath, fileName=fileName,
                                     fullOutput=fullOutput, printState=False, useGraph=use_graph,
                                     fastRelaxation=fast_relaxation, currentMap=current_map,
                                     dbg=dbg, plotCurrentMaps=plot_current_map, plotBinaryCurrentMaps=plot_binary_current_map, resume=resume,
                                     checkSteadyState=False, superconducting=sc, gap=gap, leaping=leaping,
                                     modifyR=modifyR, constQ=constQ)
    saveParameters(savePath, fileName, options, array_params)

    if dbg:
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())

    exit(0)
