import numpy as np
import matplotlib
import glob
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size' : 20}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)
from optparse import OptionParser
from matplotlib import colors
from scipy import integrate
from ast import literal_eval
import os, sys
import re
from sklearn.mixture import GaussianMixture
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from mpl_toolkits.axes_grid1 import make_axes_locatable


EPS = 1e-6
FIGSIZE=(30,16)
SCORE_NAMES=['hysteresisArea', 'jumpSeparationUp', 'jumpSeparationDown', 'jumpSeparationUpDown', 'jumpHeightUp',
             'jumpHeightDown', 'thresholdVoltageUp', 'thresholdVoltageDown', 'jumpNumUp', 'jumpNumDown', 'firstJumpV',
             'smallVPowerUp', 'smallVPowerDown', 'smallVCoefficientUp', 'smallVCoefficientDown', 'maxJumpHeight',
             'pathNum']
ACTIONS = ['plot_results', 'plot_score_by_parameter', 'plot_score_by_disorder',
           'plot_score_by_average', 'plot_scores_together', 'plot_RT', 'plot_jumps', 'plot_current_by_path',
           'IV_different_temperatures', 'perpendicular_current', 'plot_small_v_fit']
SCORE_VAL = 'score'
SCORE_HIGH_ERR = 'high_err'
SCORE_LOW_ERROR = 'low_err'


def add_text_upper_left_corner(ax, text):
    x_min,x_max,y_min,y_max = ax.axis()
    x = x_min + (x_max-x_min)/100
    y = y_max - (y_max-y_min)/9
    ax.text(x, y, text, fontsize=30)
    return ax

def add_subplot_axes(fig, ax,rect):
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]
    subax = fig.add_axes([x,y,width,height])
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax


def flattenToColumn(a):
    return a.reshape((a.size, 1))

def flattenToRow(a):
    return a.reshape((1, a.size))

class MissingFilesException(Exception):
    pass


class SingleResultsProcessor:
    """ Used for proccessing result from one dot array simulation"""
    def __init__(self, directory, file_name, fullOutput=False, vertCurrent=False, reAnalyze=True, graph=False, IT=False):
        self.basePath = os.path.join(directory, file_name)
        self.fileName = file_name
        self.directory = directory
        self.full = fullOutput
        self.vert = vertCurrent
        self.I = None
        self.V = None
        self.T = None
        self.IErr = None
        self.Imaps = None
        self.n = None
        self.Q = None
        self.nErr = None
        self.QErr = None
        self.full_I = None
        self.alternativeV = None
        self.alternativeI = None
        self.jumpScores = None
        self.smallVScores = None
        self.path_dict = None
        self.mid_idx = 0
        self.runningParams = dict()
        self.arrayParams = dict()
        self.graph = graph
        self.symmetricV = False
        try:
            self.load_results(reAnalyze, IT)
        except FileNotFoundError as e:
            raise MissingFilesException
        if not graph:
            self.load_params()
            self.rows = self.runningParams["M"]
            self.columns = self.runningParams["N"]
            self.normalize()
            if "vSym" in self.runningParams and self.runningParams["vSym"]:
                self.symmetricV = True

    def normalize(self):
        Vnorm = self.calc_param_average("C")
        Inorm = Vnorm*self.calc_param_average("R")
        if self.V is not None:
           self.V *= Vnorm
        self.I *= Inorm
        if self.IErr is not None:
            self.IErr *= Inorm
        if self.full_I is not None:
            self.full_I *= Inorm
        if self.alternativeI is not None:
            self.alternativeI *= Inorm
        if self.alternativeV is not None:
            self.alternativeV *= Vnorm
        if self.T is not None:
            self.T *= Vnorm
        if self.Imaps is not None:
            self.Imaps *= Inorm

    def reAnalyzeI(self):
        full_I = self.full_I.T if self.full_I is not None else self.I
        self.alternativeI = []
        self.alternativeIErr = []
        self.alternativeV = []
        for i in range(len(full_I)):
            min_idx = max(0, i-10)
            max_idx = min(len(full_I)-1, i+10)
            point = np.hstack([full_I[j] for j in range(min_idx, max_idx)])
            bins = bistabilityAnalysis(point)
            # bins = [point]
            if i > 0 and len(bins) > 1 and len(bins[0]) > 2 and len(bins[1]) > 2:
                avg1 = np.average(bins[0])
                avg2 = np.average(bins[1])
                # if np.abs(avg1 - self.I[i-1]) < np.abs(avg2 - self.I[i-1]):
                if len(bins[0]) > len(bins[1]):
                    self.I[i] = avg1
                    self.IErr[i] = np.std(bins[0])
                    self.alternativeI.append(avg2)
                    self.alternativeIErr.append(np.std(bins[1]))
                    self.alternativeV.append(self.V[i])
                else:
                    self.I[i] = avg2
                    self.IErr[i] = np.std(bins[1])
                    self.alternativeI.append(avg1)
                    self.alternativeIErr.append(np.std(bins[0]))
                    self.alternativeV.append(self.V[i])
            else:
                if len(bins) > 1:
                    if len(bins[0]) > len(bins[1]):
                        bin = bins[0]
                        alternativeBin = bins[1]
                    else:
                        bin = bins[1]
                        alternativeBin = bins[0]
                    if len(alternativeBin) > 0:
                        self.alternativeI.append(np.average(alternativeBin))
                        self.alternativeIErr.append(np.std(alternativeBin))
                        self.alternativeV.append(self.V[i])
                else:
                    bin = bins[0]
                self.I[i] = np.average(bin)
                self.IErr[i] = np.std(bin)
            if i == self.mid_idx:
                self.alternative_mid_idx = len(self.alternativeI)
        self.alternativeI = np.array(self.alternativeI)
        self.alternativeV = np.array(self.alternativeV)
        self.alternativeIErr = np.array(self.alternativeIErr)
        if self.Imaps is not None:
            kernel = np.ones((1000,))
            self.Imaps = np.apply_along_axis(lambda m: np.convolve(m, kernel, mode='same'), axis=0, arr=self.Imaps)
        return True


    def load_results(self, reAnalyze=False, IT=False):
        I_file = os.path.join(self.basePath + '_I.npy')
        if IT:
            T_file = os.path.join(self.basePath + '_T.npy')
            self.T = np.load(T_file)
        else:
            V_file = os.path.join(self.basePath + '_V.npy')
            self.V = np.load(V_file)
        IErr_file = os.path.join(self.basePath + '_IErr.npy')
        self.I = np.load(I_file)
        if self.graph or not os.path.isfile(IErr_file):
            self.IErr = np.zeros(self.I.shape)
        else:
            self.IErr = np.load(IErr_file)
        self.mid_idx = self.T.size //2 if IT else self.V.size // 2
        if self.full:
            n_file = os.path.join(self.basePath + '_n.npy')
            Q_file = os.path.join(self.basePath + '_Q.npy')
            nErr_file = os.path.join(self.basePath + '_nErr.npy')
            QErr_file = os.path.join(self.basePath + '_QErr.npy')
            full_I_file = os.path.join(self.basePath + '_full_I.npy')
            self.n = np.load(n_file)
            self.Q = np.load(Q_file)
            if self.graph or not os.path.isfile(QErr_file) or not os.path.isfile(nErr_file):
                self.nErr = np.zeros(self.n.shape)
                self.QErr = np.zeros(self.Q.shape)
            else:
                self.nErr = np.load(nErr_file)
                self.QErr = np.load(QErr_file)
                self.full_I = np.load(full_I_file)
            if reAnalyze:
                self.reAnalyzeI()
        if self.vert:
            vertI_file = os.path.join(self.basePath + '_vertI.npy')
            vertIErr_file = os.path.join(self.basePath + '_vertIErr.npy')
            self.vertI = np.load(vertI_file)
            self.vertIErr = np.load(vertIErr_file)
        Imaps_file = os.path.join(self.basePath + '_Imap.npy')
        if os.path.isfile(Imaps_file):
            self.Imaps = np.load(Imaps_file)
        return True

    def save_re_analysis(self):
        np.save(self.basePath + '_I_clean', self.I)
        np.save(self.basePath + '_I_Err_clean', self.IErr)
        fig = plt.figure(figsize=FIGSIZE)
        plt.errorbar(self.V[:self.mid_idx], self.I[:self.mid_idx], fmt='r.', yerr=self.IErr[:self.mid_idx],
                     errorevery=10, label="increasing voltage")
        plt.errorbar(self.V[self.mid_idx:], self.I[self.mid_idx:], fmt='b.', yerr=self.IErr[self.mid_idx:],
                     errorevery=10, label="decreasing voltage")
        if (self.I > 0.0001).any():
            plt.xlim(np.min(self.V[self.I > 0.0001]) - 0.01, np.max(self.V))
        plt.xlabel('Voltage')
        plt.ylabel('Current')
        plt.legend()
        plt.savefig(self.basePath + '_IV_clean.png')
        plt.close(fig)

    def load_params(self):
        params_file = os.path.join(self.directory,"runningParameters_" + self.fileName + ".txt")
        if not os.path.isfile(params_file):
            raise MissingFilesException
        with open(params_file, 'r') as f:
            start = False
            for line in f:
                if ' = ' in line:
                    splitted = line.split(' = ')
                    key = splitted[0]
                    try:
                        value = literal_eval(splitted[1].rstrip('\n'))
                    except Exception:
                        value = splitted[1].rstrip('\n')
                    self.runningParams[key] = value
                elif ': ' in line:
                    if start:
                        key = splitted[0]
                        try:
                            splitted[1] = splitted[1].rstrip('\n')
                            splitted[1] = re.sub('\[\s+', '[', splitted[1])
                            splitted[1] = re.sub('\s+\]', ']', splitted[1])
                            splitted[1] = re.sub('\s+', ',', splitted[1])
                            value = literal_eval(splitted[1])
                        except Exception:
                            value = splitted[1].rstrip('\n')
                        if key == "Cv":
                            value = value[1:]
                        self.arrayParams[key] = value
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
            if key == "Cv":
                value = value[1:]
            self.arrayParams[key] = value
        return True

    def get_array_param(self, key):
        return self.arrayParams[key]

    def get_running_param(self, key):
        if key in self.runningParams:
            return self.runningParams[key]
        else:
            return 0

    def getNprime(self, VL, VR):
        Ch = np.array(self.get_array_param("Ch"))
        left_part_n_prime = np.copy(Ch[:, :-1])
        left_part_n_prime[:, 1:] = 0
        left_part_n_prime = flattenToRow(left_part_n_prime)
        right_part_n_prime = np.copy(Ch[:, 1:])
        right_part_n_prime[:, :-1] = 0
        right_part_n_prime = flattenToRow(right_part_n_prime)

        return self.n + (flattenToColumn(VL).dot(left_part_n_prime) + flattenToColumn(VR).dot(right_part_n_prime)).reshape(self.n.shape)

    def calc_hysteresis_score(self, I, IErr, V, mid_idx):
        IplusErr = np.copy(I)
        IplusErr[:mid_idx] -= IErr[:mid_idx]
        IplusErr[mid_idx:] += IErr[mid_idx:]
        IminusErr = np.copy(I)
        IminusErr[:mid_idx] += IErr[:mid_idx]
        IminusErr[mid_idx:] -= IErr[mid_idx:]
        score = -integrate.trapz(I, V)
        low_err = -integrate.trapz(IplusErr, V)
        high_err = -integrate.trapz(IminusErr, V)
        return score, high_err-score, score-low_err

    def calc_threshold_voltage_up(self, I, IErr, V, mid_idx):
        I = np.array(I[:mid_idx])
        IErr = np.clip(IErr[:mid_idx],a_min=0.01, a_max=np.inf)
        matchingV = np.array(V[:mid_idx])
        dV = np.abs(self.V[1] - self.V[0])
        if (I < 2*IErr).any():
            small = np.max(matchingV[I < 2*IErr])
            big = np.max(matchingV[I < 4*IErr])
            avg = (small + big)/2
            return avg, max(big - avg, dV), max(avg - small, dV)
        else:
            return 0,0,0

    def calc_threshold_voltage_down(self, I, IErr, V, mid_idx):
        I = np.array(I[mid_idx:])
        IErr = np.clip(IErr[mid_idx:],a_min=0.01, a_max=np.inf)
        matchingV = np.array(V[mid_idx:])
        dV = np.abs(self.V[1] - self.V[0])
        if (I < 2 * IErr).any():
            small = np.max(matchingV[I < 2 * IErr])
            big = np.max(matchingV[I < 4 * IErr])
            avg = (small + big) / 2
            return avg, max(big - avg, dV), max(avg - small, dV)
        else:
            return 0, 0, 0
    def calc_number_of_jumps(self):
        _,_,_,_, Vjumps1, Vjumps2,_,_,_,_ = self.calc_jumps_freq(self.I, self.IErr)
        return len(Vjumps1) + len(Vjumps2)

    def plot_power(self):
        power = self.I*self.V
        powerPlusErr = (self.I + self.IErr)*self.V
        powerMinusErr = (self.I - self.IErr)*self.V
        plt.figure()
        plt.plot(self.V[:self.mid_idx], power[:self.mid_idx], 'b.',
                 self.V[self.mid_idx:], power[self.mid_idx:], 'r.',
                 self.V[:self.mid_idx], powerPlusErr[:self.mid_idx], 'b--',
                 self.V[self.mid_idx:], powerPlusErr[self.mid_idx:], 'r--',
                 self.V[:self.mid_idx], powerMinusErr[:self.mid_idx], 'b--',
                 self.V[self.mid_idx:], powerMinusErr[self.mid_idx:], 'r--',
                 label="horizontal")

        if self.vert:
            vertV = self.runningParams["VU"] - self.runningParams["VD"]
            vert_power = self.vertI * vertV
            vert_powerPlusErr = (self.vertI + self.vertIErr) * vertV
            vert_powerMinusErr = (self.vertI - self.vertIErr) * vertV
            plt.plot(self.V[:self.mid_idx], vert_power[:self.mid_idx], 'b.',
                     self.V[self.mid_idx:], vert_power[self.mid_idx:], 'r.',
                     self.V[:self.mid_idx], vert_powerPlusErr[:self.mid_idx], 'b--',
                     self.V[self.mid_idx:], vert_powerPlusErr[self.mid_idx:], 'r--',
                     self.V[:self.mid_idx], vert_powerMinusErr[:self.mid_idx], 'b--',
                     self.V[self.mid_idx:], vert_powerMinusErr[self.mid_idx:], 'r--',
                     label="vertical")
            plt.legend()
        plt.xlabel('Voltage')
        plt.ylabel('Power')

    def plot_resistance(self, out=None, show=False):
        Vup = self.V[:self.mid_idx]
        Iup = self.I[:self.mid_idx]
        Vup = Vup[Iup != 0]
        Iup = Iup[Iup != 0]
        resistance_up = Vup/Iup

        Vdown = self.V[self.mid_idx:]
        Idown = self.I[self.mid_idx:]
        Vdown = Vdown[Idown != 0]
        Idown = Idown[Idown != 0]
        resistance_down = Vdown / Idown

        fig = plt.figure(figsize=FIGSIZE)
        plt.semilogy(Vup, resistance_up, 'r.',label="horizontal resistance increasing v")
        plt.semilogy(Vdown, resistance_down, 'b.', label="horizontal resistance decreasing v")

        if self.vert:
            vertV = 0.2
            vertIup = self.vertI[:self.mid_idx]
            Vup = self.V[:self.mid_idx]
            Vup = Vup[vertIup != 0]
            vertIup = vertIup[vertIup != 0]
            resistance_up = vertV/vertIup
            vertIdown = self.vertI[self.mid_idx:]
            Vdown = self.V[self.mid_idx:]
            Vdown = Vdown[vertIdown != 0]
            vertIdown = vertIdown[vertIdown != 0]
            resistance_down = vertV / vertIdown

            plt.semilogy(Vup, resistance_up, 'm.', label="vertical resistance increasing v")
            plt.semilogy(Vdown, resistance_down, 'g.', label="vertical resistance dencreasing v")
            plt.legend()
        plt.xlabel('Voltage')
        plt.ylabel('Resistance')
        if out is not None:
            plt.savefig(out + "_log_resistance.png")
        if show:
            plt.show()
        else:
            plt.close(fig)

    def plot_jumps_freq(self, by_occupation=False, by_path_occupation=False, ax1=None, fig=None):
        V = self.V[self.V>=0]
        I = self.I[self.V>=0]
        IErr = self.IErr[self.V>=0]
        mid_idx = np.argmax(V)
        thresholdV = np.min(V[I>EPS])
        IV_diffV1, IV_diffV2, IV_diff1, IV_diff2, IV_Vjumps1, IV_Vjumps2, IV_freq1, IV_freq2, IV_fou1, IV_fou2 =\
            self.calc_jumps_freq(I[V >= thresholdV], IErr[V>=thresholdV], V=V[V>=thresholdV], mid_idx=np.argmax(V[V>=thresholdV]),
                                 threshold_factor=2)
        if ax1 is None:
            fig, ax1 = plt.subplots(figsize=FIGSIZE)
        fig2 = None
        fig3 = None
        island_colors = None
        ax1.set_xlabel('$V\\frac{\\left<C\\right>}{e}$')
        ax1.set_ylabel('$I\\frac{\\left<C\\right>\\left<R\\right>}{e}$')
        ax1.plot(V[:mid_idx], I[:mid_idx], 'r*', zorder=1)
        ax1.plot(V[mid_idx:], I[mid_idx:], 'b*', zorder=2)
        # if not self.graph:
            # ax1.set_title("$\\frac{\\left<C_G\\right>}{\\left<C\\right>} = " +
            #               str(np.round(self.calc_param_average("CG")/self.calc_param_average("C")*100)/100) + "$",
            #               position=(0.5,0.9))
        arrow_length = np.max(I)/60
        arrow_head_width = np.max(V)/40
        if self.Imaps is not None and self.path_dict is None:
            self.path_dict = self.calc_paths()
            colors = matplotlib.cm.gist_rainbow(np.linspace(0, 1, len(self.path_dict.keys())))
        if (self.n is None and self.Imaps is None) or not self.full:
            for v, diff in zip(IV_Vjumps1, IV_diff1):
                i = I[V == v][0]
                ax1.arrow(v, i, 0, diff, fc='r',ec='r', head_length=arrow_head_width, head_width=arrow_head_width, zorder=3)
            for v, diff in zip(IV_Vjumps2, IV_diff2):
                i = I[V == v][-1]
                ax1.arrow(v, i, 0, -diff,fc='b',ec='b', head_length=arrow_head_width, head_width=arrow_head_width, zorder=3)
        else:
            ax2 = add_subplot_axes(fig, ax1, [0.7,0.1,0.3,0.4])
            fig3, fou_ax = plt.subplots(1, 2, figsize=FIGSIZE)
            for ax in fou_ax:
                ax.set_xlabel("Frequency [$\\frac{e}{V\\left<C\\right>}$]")
                ax.set_ylabel("Fourier Amplitude [a.u.]")
            fou_ax[0].set_title("Increasing voltage")
            fou_ax[1].set_title("Decreasing voltage")
            fou_ax[0].plot(IV_freq1, IV_fou1)
            fou_ax[1].plot(IV_freq2, IV_fou2)
            if not by_occupation and self.Imaps is not None:
                ax2.set_ylabel("$I\\frac{\\left<C\\right>\\left<R\\right>}{e}$", fontsize=20)
                ax2.yaxis.set_label_position("right")
                ax2.yaxis.tick_right()
                if len(self.path_dict.keys()) > 0:
                    fig2, hist_axes = plt.subplots(4, len(self.path_dict.keys())+1, figsize=FIGSIZE)
                    hist_axes[0][0].hist(IV_diff1)
                    hist_axes[1][0].hist(IV_diff2)
                    hist_axes[2][0].hist(IV_diffV1)
                    hist_axes[3][0].hist(IV_diffV2)

                for index,path in enumerate(sorted(self.path_dict.keys())):
                    x = np.array(self.path_dict[path][1])
                    if len(x) < 10:
                        continue
                    xerr = np.zeros(x.shape)
                    V_path = np.array(self.path_dict[path][0])
                    mid_idx_path = np.argmax(V_path)
                    color = colors[index]
                    diffV1, diffV2, diff1, diff2, Vjumps1, Vjumps2, freq1, freq2, fou1, fou2 =\
                    self.calc_jumps_freq(x, xerr, V=V_path, mid_idx=mid_idx_path,threshold_factor=10,
                                         window_size_factor=1, absolute_threshold=0.003)
                    hist_axes[0][index+1].hist(diff1, color=color)
                    hist_axes[1][index+1].hist(diff2, color=color)
                    hist_axes[2][index+1].hist(diffV1, color=color)
                    hist_axes[3][index+1].hist(diffV2, color=color)
                    fou_ax[0].plot(freq1, fou1, color=color)
                    fou_ax[1].plot(freq2, fou2, color=color)
                    ax2.plot(V_path[:mid_idx_path], x[:mid_idx_path], color=color, linestyle='dashed',zorder=1)
                    v_open_index, = np.where(V == V_path[0])
                    if v_open_index[0] > 0:
                        ax1.arrow(V[v_open_index[0]-1], I[v_open_index[0]-1], 0,x[0],
                                  fc=color, ec=color, head_length=2*arrow_length, head_width=2*arrow_head_width,zorder=3)
                    ax1.arrow(V_path[-1], I[V == V_path[-1]][-1], 0, -x[-1], fc=color, ec=color,
                              head_length=2 * arrow_length,
                              head_width=2*arrow_head_width, zorder=3)
                    for v, diff in zip(Vjumps1, diff1):
                        i = I[V == v][0]
                        ax1.arrow(v, i, 0, diff, fc=color, ec=color, head_length=arrow_length,
                                  head_width=arrow_head_width, zorder=3)
                    ax2.plot(V_path[mid_idx_path:], x[mid_idx_path:], color=color, linestyle='dotted',zorder=2)
                    for v, diff in zip(Vjumps2, diff2):
                        i = I[V == v][-1]
                        ax1.arrow(v, i, 0, -diff, fc=color, ec=color, head_length=arrow_length,
                                  head_width=arrow_head_width, zorder=3)

            if self.n is not None and (by_occupation or by_path_occupation):
                VL = self.V / 2 if self.symmetricV else self.V
                VR = -self.V / 2 if self.symmetricV else np.zeros(self.V.shape)
                n = self.getNprime(VL, VR).reshape((self.n.shape[0], self.n.shape[1] * self.n.shape[2]))
                n = n[self.V >= 0]
                nErr = self.nErr[self.V >= 0].reshape(n.shape)
                if by_path_occupation and self.path_dict is not None:
                    current_ax = add_subplot_axes(fig, ax1, [0, 0.12, 0.25, 0.3])
                    current_ax.yaxis.set_label_position("right")
                    current_ax.yaxis.tick_right()
                    current_ax.set_ylabel("$\\sum\\left<n'\\right>$",fontsize=20)
                    path_n = np.zeros((n.shape[0],len(self.path_dict)))
                    for index,path in enumerate(sorted(self.path_dict.keys())):
                        for i,j in path:
                            if j < self.columns:
                                path_n[:,index] += n[:,i*self.columns + j]
                    index_range = len(self.path_dict)
                    n = path_n
                else:
                    index_range = self.columns*self.rows
                    current_ax = ax2
                    current_ax.set_ylabel("$\\left<n'\\right>$")
                colors = matplotlib.cm.gist_rainbow(np.linspace(0, 1, index_range))
                for index in range(index_range):
                    color = colors[index]
                    x = n[:,index]
                    xerr = nErr[:,index]
                    diffV1, diffV2, diff1, diff2, Vjumps1, Vjumps2, freq1, freq2, fou1, fou2 = \
                        self.calc_jumps_freq(x, xerr, V=V,mid_idx=mid_idx, up_and_down=True, threshold_factor=10,
                                             window_size_factor=1, absolute_threshold=0.01)
                    current_ax.plot(V[:mid_idx], x[:mid_idx], color=color, linestyle='dashed')
                    current_ax.plot(V[mid_idx:], x[mid_idx:], color=color, linestyle='dotted')
                    fou1 *= np.max(IV_fou1)/np.max(fou1)
                    fou2 *= np.max(IV_fou2)/np.max(fou2)
                    fou_ax[0].plot(freq1, fou1, color=color)
                    fou_ax[1].plot(freq2, fou2, color=color)
                    if by_occupation:
                        for v, diff in zip(Vjumps1, diff1):
                            i = I[V==v][0]
                            ax1.arrow(v, i, 0, diff/50, fc=color, ec=color, head_length=arrow_length,
                                      head_width=arrow_head_width)
                        for v, diff in zip(Vjumps2, diff2):
                            i = I[V==v][-1]
                            ax1.arrow(v, i, 0, -diff/50, fc=color, ec=color,
                                      head_length=arrow_length, head_width=arrow_head_width)
                        island_colors = colors
            ax3_location = [0.05,0.6,0.5,0.2] if self.rows < 3 else [0.05,0.5,0.6,0.5]
            ax3 = add_subplot_axes(fig, ax1, ax3_location)
            self.plot_array_params("RC", ax3, fig, island_colors=island_colors)
            if self.Imaps is not None and not by_occupation:
                plot_paths(sorted(self.path_dict.keys()), ax3, colors)
        return fig, fig2, fig3


    def calc_jumps_freq(self, x, xerr, V=None, mid_idx=None, threshold_factor=1, window_size_factor=1,
                        absolute_threshold=0.001, up_and_down=False):
        if V is None:
            V = self.V
        if mid_idx is None:
            mid_idx = self.mid_idx
        vstep = self.get_running_param("vStep")
        if vstep == 0:
            vstep = 0.01
        window_size = int(np.round(0.05*window_size_factor / vstep))
        diff1 = average_diff(x[:mid_idx], window_size)
        diff2 = average_diff(x[mid_idx:], window_size)
        if up_and_down:
            diff1 = np.abs(diff1)
            diff2 = np.abs(diff2)
        dV = self.V[1]-self.V[0]
        if diff1.size > 0:
            freq1, fou1 = self.calc_fourier(dV, diff1)
        else:
            freq1, fou1 = [],[]
        if diff2.size > 0:
            freq2, fou2 = self.calc_fourier(dV, diff2)
        else:
            freq2, fou2 = [], []
        if np.sum(diff1) < 0:
            diff1 *= -1
        if np.sum(diff2) < 0:
            diff2 *= -1

        if len(diff1) < 1:
            threshold1 = 0.1
        else:
            small_diff1 = diff1[np.logical_and(diff1 > 0, diff1 < np.average(diff1))]
            if len(small_diff1) < 1:
                threshold1 = 0.1
            else:
                threshold1 = max(threshold_factor * np.std(small_diff1),absolute_threshold)
        if len(diff2) < 1:
            threshold2 = 0.1
        else:
            small_diff2 = diff2[np.logical_and(diff2 > 0, diff2 < np.average(diff2))]
            if len(small_diff2) < 1:
                threshold2 = 0.1
            else:
                threshold2 = max(threshold_factor * np.std(small_diff2),absolute_threshold)
        jumps_indices1, _ = find_peaks(diff1, prominence=0.5*threshold1, height=threshold1)
        height1 = diff1[jumps_indices1]
        Vup = V[:mid_idx]
        Vjumps1 = Vup[jumps_indices1]
        jumps_indices2, _ = find_peaks(diff2, prominence=0.5*threshold2, height=threshold2)
        height2 = diff2[jumps_indices2]
        Vdown = V[mid_idx:]
        Vjumps2 = Vdown[jumps_indices2]
        diffV1 = np.diff(Vjumps1)
        diffV2 = np.diff(Vjumps2)
        return diffV1, diffV2, height1, height2, Vjumps1, Vjumps2, freq1, freq2, fou1, fou2

    def get_Vjumps(self):
        _, _, _, _, Vjumps1, Vjumps2,_,_,_,_ = self.calc_jumps_freq(self.I, self.IErr)
        return Vjumps1, Vjumps2

    def calc_jumps_scores(self, I, IErr, V, mid_idx):
        diffV1, diffV2, diff1, diff2, Vjumps1, Vjumps2,_,_,_,_ = self.calc_jumps_freq(I, IErr, V=V, mid_idx=mid_idx,
                                                                                      window_size_factor=1,
                                                                                      threshold_factor=1,
                                                                                      absolute_threshold=0.001)
        self.jumpScores = dict()
        self.jumpScores['jumpSeparationUp'] = (np.average(diffV1), np.std(diffV1), np.std(diffV1)) if diffV1.size > 0 else (0,0,0)
        self.jumpScores['jumpSeparationDown'] = (-np.average(diffV2), -np.std(diffV2), -np.std(diffV2)) if diffV2.size > 0 else (0,0,0)
        self.jumpScores['jumpHeightUp'] = (np.average(diff1), np.std(diff1), np.std(diff1)) if diff1.size > 0 else (0,0,0)
        self.jumpScores['maxJumpHeight'] = (diff1[0], 0,0) if diff1.size > 0 else (0,0,0)
        self.jumpScores['jumpHeightDown'] = (np.average(diff2), np.std(diff2), np.std(diff2)) if diff2.size > 0 else (0,0,0)
        self.jumpScores['jumpNumUp'] = (Vjumps1.size, 0, 0)
        self.jumpScores['jumpNumDown'] = (Vjumps2.size, 0, 0)
        if len(Vjumps2) > 0 and len(Vjumps1)> 0:
            down = Vjumps2[-1]
            up_jumps_larger_than_down = Vjumps1[Vjumps1>down]
            if len(up_jumps_larger_than_down) > 0:
                up = np.min(up_jumps_larger_than_down)

            else:
                up = down
            self.jumpScores['jumpSeparationUpDown'] = ((up-down), 0, 0)
            self.jumpScores['firstJumpV'] = (up,0,0)
        else:
            self.jumpScores['jumpSeparationUpDown'] = (0, 0, 0)
            self.jumpScores['firstJumpV'] = (0,0,0)

    def calc_small_v_scores(self, I, IErr, V, mid_idx):
        self.smallVScores = dict()
        paramsUp, covUp, paramsDown, covDown = self.power_law_fit(I, IErr, V, mid_idx)
        self.smallVScores['smallVPowerUp'] = paramsUp[1], covUp[1], covUp[1]
        self.smallVScores['smallVPowerDown'] = paramsDown[1], covDown[1], covDown[1]
        self.smallVScores['smallVCoefficientUp'] = paramsUp[0], covUp[0], covUp[0]
        self.smallVScores['smallVCoefficientDown'] = paramsDown[0], covDown[0], covDown[0]


    def calc_fourier(self, dx, y):
        fou = np.abs(np.fft.rfft(y))
        freq = np.fft.rfftfreq(y.size, dx)
        fou[0] = 0
        return freq, fou

    def plot_diff_conductance(self):
        dI1 = np.diff(self.I[:self.mid_idx])
        dI2 = np.diff(self.I[self.mid_idx:])
        dV1 = np.diff(self.V[:self.mid_idx])
        dV2 = np.diff(self.V[self.mid_idx:])
        plt.figure()
        plt.plot(self.V[:self.mid_idx-1], dI1/dV1, label="increasing voltage")
        plt.plot(self.V[self.mid_idx:-1], dI2/dV2, label="decreasing voltage")
        plt.xlabel("Voltage")
        plt.ylabel("Conductivity")
        plt.legend()

    def plot_conductance(self):
        condUp = self.I[1:self.mid_idx]/self.V[1:self.mid_idx]
        condDown = self.I[self.mid_idx:-1]/self.V[self.mid_idx:-1]
        plt.figure()
        plt.plot(self.V[1:self.mid_idx], condUp, label="increasing voltage")
        plt.plot(self.V[self.mid_idx:-1], condDown, label="decreasing voltage")
        if self.vert:
            vertV = self.runningParams["VU"] - self.runningParams["VD"]
            vertcondUp = np.abs(self.vertI[:self.mid_idx]) / vertV
            vertcondDown = np.abs(self.vertI[self.mid_idx:]) / vertV
            plt.plot(self.V[:self.mid_idx], vertcondUp, label="vertical, increasing voltage")
            plt.plot(self.V[self.mid_idx:], vertcondDown, label="vertical, decreasing voltage")

        plt.xlabel("Voltage")
        plt.ylabel("Conductivity")
        plt.legend()

    def calc_resistance_hysteresis_score(self, line):
        R = np.array(self.arrayParams["Rh"])
        C = np.array(self.arrayParams["Ch"])
        Cv = np.array(self.arrayParams["Cv"])
        # chosen_line = 0
        # max_min_c = 0
        # for line in range(1,C.shape[0]-1):
        #     c = C[line,:-1] + C[line,1:] + Cv[line, :] + Cv[line-1,:]
        #     min_c = np.min(c)
        #     if min_c > max_min_c:
        #         max_min_c = min_c
        #         chosen_line = line
        min_resistance = np.inf
        for line in range(1, R.shape[0]-1):
            resistance = np.sum(R[line,:])
            if resistance < min_resistance:
                min_resistance = resistance
                chosen_line = line
        line = chosen_line
        R = R[line, :]
        C = C[line, :]
        Cup = Cv[line, :]
        Cdown = Cv[line+1,:]
        score = 0
        mid = int(R.size/2)
        for i in range(mid):
            exit = R[i+1]
            enter = np.max(R[:i+1])
            # score += (mid - i)*(exit/enter)/(C[i] + C[i+1])
            score += (exit /enter) / (C[i] + C[i + 1] + Cup[i] + Cdown[i])**2
        for i in range(mid, R.size-1):
            exit = np.max(R[i + 1:])
            enter = R[i]
            # score += (1+i - mid)*(enter / exit)/(C[i] + C[i+1])
            score += (enter/exit) / (C[i] + C[i + 1]+ Cup[i] + Cdown[i])**2
        return score

    def calc_score(self, scoreName, I=None, IErr=None, V=None, mid_idx=None, re_calculate=False):
        if I is None:
            I = self.I
        if IErr is None:
            IErr = self.IErr
        if V is None:
            V = self.V
        if mid_idx is None:
             mid_idx = self.mid_idx
        if 'jump' in scoreName or 'Jump' in scoreName:
            if self.jumpScores is None or re_calculate:
                self.calc_jumps_scores(I, IErr, V, mid_idx)
            return self.jumpScores[scoreName]
        if 'smallV' in scoreName:
            if self.smallVScores is None or re_calculate:
                self.calc_small_v_scores(I, IErr, V, mid_idx)
            return self.smallVScores[scoreName]
        elif scoreName == "hysteresisArea":
            return self.calc_hysteresis_score(I, IErr, V, mid_idx)
        elif scoreName == 'thresholdVoltageUp':
            return self.calc_threshold_voltage_up(I, IErr, V, mid_idx)
        elif scoreName == 'thresholdVoltageDown':
            return self.calc_threshold_voltage_down(I, IErr, V, mid_idx)
        elif scoreName == 'pathNum':
            return self.calc_path_num()
        else:
            raise NotImplementedError

    def calc_path_num(self):
        if self.Imaps is None:
            return 0,0,0
        if self.path_dict is None:
            self.path_dict = self.calc_paths()
        return len(self.path_dict.keys()),0,0

    def calc_score_for_paths(self, score_name):
        if self.path_dict is None:
            self.path_dict = self.calc_paths()
        path_scores = dict()
        for path in self.path_dict:
            V = np.array(self.path_dict[path][0])
            I = np.array(self.path_dict[path][1])
            IErr = np.zeros((len(I),))
            mid_idx = np.argmax(V)
            path_scores[path] = self.calc_score(score_name,I=I, IErr=IErr, V=V, mid_idx=mid_idx, re_calculate=True)
        return path_scores

    def print_scores(self, for_paths=False):
        print("Scores for file %s, from dir %s" % (self.fileName, self.directory))
        for score in SCORE_NAMES:
            score_val, high, low = self.calc_score(score)
            print("%s: %f (%f, %f)" % (score, score_val, score_val - low, score_val + high))
        if for_paths:
            for score in SCORE_NAMES:
                path_scores = self.calc_score_for_paths(score)
                for path in path_scores:
                    score_val, low, high = path_scores[path]
                    print("For path %s score %s is %f (%f, %f)" % (path, score, score_val,
                                                                   score_val - low, score_val + high))

    def plot_array_params(self, parameter,ax=None, fig=None, island_colors=None):
        M = self.runningParams["M"]
        N = self.runningParams["N"]
        connection = np.zeros((M*6-4,N*6+4))
        C_avg = self.calc_param_average("C")
        R_avg = self.calc_param_average("R")
        if parameter in ["R", "C"]:
            horzRows = np.repeat(np.arange(0, M * 6 - 4, 6), 2)
            horzRows[1::2] +=1
            horzCols = np.repeat(np.arange(0, 6 * N + 4, 6), 4)
            horzCols[1::4] += 1
            horzCols[2::4] += 2
            horzCols[3::4] += 3
            vertRows = np.repeat(np.arange(2, M * 6 -4, 6), 4)
            vertRows[1::4] += 1
            vertRows[2::4] += 2
            vertRows[3::4] += 3
            vertCols = np.repeat(np.arange(4, 6 * N + 4, 6), 2)
            vertCols[1::2] += 1
            if parameter == "C":
                data_horz = np.array(self.arrayParams["Ch"])/C_avg
                data_vert = np.array(self.arrayParams["Cv"])/C_avg
            elif parameter == "R":
                data_horz = np.array(self.arrayParams["Rh"]) / R_avg
                data_vert = np.array(self.arrayParams["Rv"]) / R_avg
            connection_mask = np.ones(connection.shape)
            connection_mask[np.ix_(horzRows, horzCols)] = 0
            connection_mask[np.ix_(vertRows, vertCols)] = 0
            if parameter == "C":
                cmap = 'Greens'
            elif parameter == "R":
                cmap = 'Reds'
            data = np.ones((self.columns, self.rows))
            connection[np.ix_(horzRows, horzCols)] = np.repeat(np.repeat(data_horz,4,axis=1),2,axis=0)
            if M > 1:
                connection[np.ix_(vertRows, vertCols)] = np.repeat(np.repeat(data_vert,4,axis=0), 2, axis=1)
        elif parameter == "RC":
            connection2 = np.zeros(connection.shape)
            horzRowsR = np.arange(0, M * 6 - 4, 6)
            horzRowsC = np.arange(1, M * 6 - 4, 6)
            horzCols = np.repeat(np.arange(0, 6 * N + 4, 6), 4)
            horzCols[1::4] += 1
            horzCols[2::4] += 2
            horzCols[3::4] += 3
            vertRows = np.repeat(np.arange(2, M * 6 -4, 6), 4)
            vertRows[1::4] += 1
            vertRows[2::4] += 2
            vertRows[3::4] += 3
            vertColsR = np.arange(4, 6 * N + 4, 6)
            vertColsC = np.arange(5, 6 * N + 4, 6)
            data_horz_C = np.array(self.arrayParams["Ch"]) / C_avg
            data_vert_C = np.array(self.arrayParams["Cv"]) / C_avg
            data_horz_R = np.array(self.arrayParams["Rh"]) / R_avg
            data_vert_R = np.array(self.arrayParams["Rv"]) / R_avg
            connection_mask = np.ones(connection.shape)
            connection_mask[np.ix_(horzRowsR, horzCols)] = 0
            connection_mask[np.ix_(vertRows, vertColsR)] = 0
            connection2_mask = np.ones(connection.shape)
            connection2_mask[np.ix_(horzRowsC, horzCols)] = 0
            connection2_mask[np.ix_(vertRows, vertColsC)] = 0
            data = np.ones((self.rows, self.columns))
            if island_colors is not None:
                data = np.array([colors.to_rgba(c) for c in island_colors])
                data = data.reshape((self.rows, self.columns,4))
            connection[np.ix_(horzRowsR, horzCols)] = np.repeat(data_horz_R, 4, axis=1)
            connection2[np.ix_(horzRowsC, horzCols)] = np.repeat(data_horz_C, 4, axis=1)
            if M > 1:
                connection[np.ix_(vertRows, vertColsR)] = np.repeat(data_vert_R, 4, axis=0)
                connection2[np.ix_(vertRows, vertColsC)] = np.repeat(data_vert_C, 4, axis=0)
            cmap = 'Reds'
            cmap2 = 'Greens'
        else:
            data = self.arrayParams["CG"]/C_avg if parameter == "CG" else\
                self.arrayParams["RG"]/R_avg if parameter == "RG" else self.arrayParams["VG"]*C_avg
            data = np.reshape(data, (self.runningParams["M"], self.runningParams["N"]))
            cmap = 'RdBu' if parameter == "VG" else 'Greens'
        rows = np.repeat(np.arange(0, M* 6 -4, 6), 2)
        rows[1::2] +=1
        cols = np.repeat(np.arange(4, 6 * N + 4, 6),2)
        cols[1::2] +=1


        if island_colors is not None:
            dots = np.zeros(connection.shape + (4,))
            dots[np.ix_(rows, cols, np.arange(4))] = np.repeat(np.repeat(data, 2, axis=0), 2, axis=1)
            dots_mask = np.ones(connection.shape + (4,))
            dots_mask[np.ix_(rows, cols, np.arange(4))] = 0
        else:
            dots = np.zeros(connection.shape)
            dots[np.ix_(rows, cols)] = np.repeat(np.repeat(data, 2, axis=0), 2, axis=1)
            dots_mask = np.ones(connection.shape)
            dots_mask[np.ix_(rows, cols)] = 0
        if ax is None:
            fig, ax = plt.subplots()
        im1 = ax.imshow(np.ma.masked_array(connection, connection_mask), vmin=np.min(connection),
                       vmax=np.max(connection), cmap=cmap, aspect='equal')
        ax.imshow(np.ma.masked_array(dots, dots_mask), vmin=np.min(dots),
                  vmax=np.max(dots), cmap='Blues', aspect='equal')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        if parameter == "RC":
            im2 = ax.imshow(np.ma.masked_array(connection2, connection2_mask), vmin=np.min(connection2),
                            vmax=np.max(connection2), cmap=cmap2, aspect='equal')
            cb = fig.colorbar(im1, cax=cax)
            cax2 = divider.append_axes("right", size="5%", pad=0.7)
            cb2 = fig.colorbar(im2, cax=cax2)
            cb.set_label("$R/\\left<R\\right>$", fontsize=15)
            cb2.set_label("$C/\\left<C\\right>$", fontsize=15)
        else:
            cb = fig.colorbar(im1, cax=cax)
            norm_parameter = "C" if parameter in ["C","CG"] else "R"
            cb.set_label("$" + parameter + "/\\left<" + norm_parameter + "\\right>$")
        y_label_list = ['row ' + str(row) for row in range(self.rows)]
        if self.columns > 3:
            x_label_list = ['col ' + str(col) for col in range(0,self.columns,2)]
            ax.set_xticks(cols[::4])
        else:
            x_label_list = ['col ' + str(col) for col in range(self.columns)]
            ax.set_xticks(cols[::2])
        ax.set_yticks(rows[::2])
        ax.set_xticklabels(x_label_list)
        ax.set_yticklabels(y_label_list)
        return ax

    def createCapacitanceMatrix(self):
        """
        Creates the inverse capacitance matrix
        """
        Ch = np.array(self.get_array_param("Ch"))
        Cv = np.array(self.get_array_param("Cv"))
        if Cv.size == 0:
            self.invC = np.array([[1/np.sum(Ch)]])
            return True
        Cv = np.pad(Cv,((1,0),(0,0)))
        diagonal = Ch[:,:-1] + Ch[:,1:] + Cv + np.roll(Cv, -1, axis=0)
        second_diagonal = np.copy(Ch[:,1:])
        second_diagonal[:,-1] = 0
        second_diagonal = second_diagonal.flatten()
        second_diagonal = second_diagonal[:-1]
        n_diagonal = np.copy(Cv[1:,:])
        C_mat = np.diagflat(diagonal) - np.diagflat(second_diagonal,k=1) - np.diagflat(second_diagonal,k=-1)\
                -np.diagflat(n_diagonal, k=self.columns) - np.diagflat(n_diagonal, k=-self.columns)
        self.invC = np.linalg.inv(C_mat)

        return True
    def plot_voltage(self):
        self.createCapacitanceMatrix()
        print(np.diagonal(self.invC))
        VL = self.V / 2 if self.symmetricV else self.V
        VR = -self.V / 2 if self.symmetricV else np.zeros(self.V.shape)
        n = self.getNprime(VL, VR).reshape((self.n.shape[0], self.n.shape[1] * self.n.shape[2]))
        Q = self.Q.reshape((self.Q.shape[0], self.Q.shape[1] * self.Q.shape[2]))
        self.q = n + Q
        self.max = len(self.q)-1
        self.i=0
        self.fig = plt.figure()
        self.X = np.arange(self.columns)
        self.Y = np.arange(self.rows)
        self.X, self.Y = np.meshgrid(self.X, self.Y)
        self.fig.canvas.mpl_connect("key_press_event", self.click)
        self.plotV(0)
        plt.show()
        plt.draw()
    def plotV(self,idx):
        v = self.invC.dot(self.q[idx])
        Z = v.reshape(self.rows, self.columns)
        plt.clf()
        ax = self.fig.gca(projection='3d')
        # ax.plot_surface(self.X, self.Y, Z)
        ax.scatter(self.X, self.Y, Z)
        ax.set_title("V=" + str(self.V[idx]))
    def click(self,event):
        if event.key == "right" and self.i < self.max:
            self.i +=1
            self.plotV(self.i)
            plt.draw()
        if event.key == "left" and self.i > 0:
            self.i -= 1
            self.plotV(self.i)
            plt.draw()
        if event.key == "." and self.i < self.max - 10:
            self.i +=10
            self.plotV(self.i)
            plt.draw()
        if event.key == "," and self.i > 10:
            self.i -= 10
            self.plotV(self.i)
            plt.draw()

    def plot_differences(self):
        # plt.figure()
        I_up = self.I[:self.mid_idx]
        I_down = self.I[self.mid_idx:]
        V = self.V[:self.mid_idx]
        xmin = np.min(V[I_up > 0.001]) - 0.1
        plt.plot(V, np.flip(I_down) - I_up, label=self.directory)
        plt.xlim(xmin, np.max(V))
        plt.ylabel("I")
        plt.xlabel("V")


    def plot_IV(self, ax=None, fig=None, err=True, alternative=False, errorevery=10,
                Vnorm=1, Inorm=1, Vlabel='Voltage', Ilabel='Current', shift=0, fmt_up='r.', fmt_down='b.'):

        if ax is None:
            fig, ax = plt.subplots()
        if err:
            ax.errorbar(self.V[:self.mid_idx]/Vnorm, self.I[:self.mid_idx]/Inorm + shift, fmt=fmt_up,
                         yerr=self.IErr[:self.mid_idx]/Inorm,
                         errorevery=errorevery)

            ax.errorbar(self.V[self.mid_idx:]/Vnorm, self.I[self.mid_idx:]/Inorm + shift, fmt=fmt_down,
                         yerr=self.IErr[self.mid_idx:]/Inorm,
                         errorevery=errorevery)
        else:
            ax.plot(self.V[:self.mid_idx]/Vnorm, self.I[:self.mid_idx]/Inorm + shift, fmt_up)
            ax.plot(self.V[self.mid_idx:]/Vnorm, self.I[self.mid_idx:]/Inorm + shift, fmt_down)
        if alternative and self.alternativeV is not None:
            ax.errorbar(self.alternativeV[:self.alternative_mid_idx]/Vnorm,
                         self.alternativeI[:self.alternative_mid_idx]/Inorm + shift,
                         fmt='m.', yerr=self.alternativeIErr[:self.alternative_mid_idx]/Inorm,
                         errorevery=errorevery)
            ax.errorbar(self.alternativeV[self.alternative_mid_idx:]/Vnorm ,
                         self.alternativeI[self.alternative_mid_idx:]/Inorm + shift,
                         fmt='c.', yerr=self.alternativeIErr[self.alternative_mid_idx:]/Inorm,
                         errorevery=errorevery)

        ax.set_xlabel(Vlabel)
        ax.set_ylabel(Ilabel)
        return ax, fig

    def plot_RT(self, ax, err=True, errorevery=10, fit=True, shift=0, color='green', label=""):
        V = self.get_running_param("Vmin") * self.calc_param_average("C")*2
        T = self.T *2
        I = self.I*4
        IErr = self.IErr*4
        IErr[IErr < EPS] = EPS
        ax.set_yscale("log", nonposy='clip')
        if err:
            ax.errorbar(T[I>0], (V/I[I>0]) * shift, color=color, marker='.',linestyle="",
                         yerr=IErr[I>0]*V/I[I>0], errorevery=errorevery)
        else:
            ax.plot(T[I>0], (V/I[I>0])* shift, color=color, marker='.',linestyle="")
        if fit:
            def arrhenius(T, a, b):
                return a*np.exp(-b/T)
            def shifted_arrhenius(T, a, b, c):
                return a * np.exp(-b / T) + c
            if I[0] == 0:
                fit_param, fit_cov = curve_fit(arrhenius, T, I / V, p0=((I[-1] - I[0]) / V, 0.1),
                                               bounds=(0, np.inf), sigma=IErr)
                G, T0 = fit_param
                G_err, T0_err = np.sqrt(np.diag(fit_cov))
                Rinf = 1 / G
                Rinf_err = Rinf * G_err
                Rinf, Rinf_err = significant_figures(Rinf, Rinf_err)
                T0,T0_err = significant_figures(T0, T0_err)
                ax.text(np.max(T)-0.37, 1.9 * shift * (V / I[-1]),
                        "$R_{\\infty}/\\left(R_1+R_2\\right) = %s \\pm %s$, $T_0\\frac{C_1+C_2}{e^2} = %s \\pm %s$" % (
                            str(Rinf), str(Rinf_err), str(T0), str(T0_err)), color=color)
                fit_result = arrhenius(T, *fit_param)
            else:
                fit_param, fit_cov = curve_fit(shifted_arrhenius, T, I/V, p0=((I[-1]-I[0])/V, 2,I[3]/V), bounds=(0,np.inf),
                                               sigma=IErr)
                G, T0, G0 = fit_param
                G_err, T0_err, G0_err = np.sqrt(np.diag(fit_cov))
                Rinf = 1 / (G + G0)
                Rinf_err = Rinf * np.sqrt(G_err ** 2 + G0_err ** 2)
                R0 = 1 / G0
                R0_err = R0 * G0_err
                Rinf, Rinf_err = significant_figures(Rinf, Rinf_err)
                T0, T0_err = significant_figures(T0, T0_err)
                R0, R0_err = significant_figures(R0, R0_err)
                ax.text(np.max(T)-0.45,1.9*shift*(V/I[-1]),"$R_{0}/\\left(R_1+R_2\\right) = %s \\pm %s$, $R_{\\infty}/\\left(R_1+R_2\\right)= %s \\pm %s$, $T_0\\frac{C_1+C_2}{e^2}  = %s \\pm %s$" % (
                str(R0), str(R0_err), str(Rinf), str(Rinf_err), str(T0), str(T0_err)), color=color)
                fit_result = shifted_arrhenius(T, *fit_param)
            ax.plot(T, (1/fit_result) * shift, color=color,label=label)
        ax.set_ylim(1, np.max(V/I[I>0] *shift)*2000)
        ax.set_xlim(0, max(ax.get_xlim()[1], np.max(T) + 0.01))
        return V


    def plot_results(self):
        plt.figure(figsize=FIGSIZE)
        self.plot_IV(self.directory, err=True, alternative=True)
        if self.full:
            VL = self.V/2 if self.symmetricV else self.V
            VR = -self.V/2 if self.symmetricV else np.zeros(self.V.shape)
            n = self.getNprime(VL, VR).reshape((self.n.shape[0], self.n.shape[1] * self.n.shape[2]))
            Q = self.Q.reshape((self.Q.shape[0], self.Q.shape[1] * self.Q.shape[2]))
            nErr = self.nErr.reshape((self.nErr.shape[0], self.nErr.shape[1] *self. nErr.shape[2]))
            QErr = self.QErr.reshape((self.QErr.shape[0], self.QErr.shape[1] * self.QErr.shape[2]))
            plt.figure(figsize=FIGSIZE)
            for i in range(len(n[0])):
                plt.errorbar(self.V[:self.mid_idx], n[:self.mid_idx, i] + i, yerr=nErr[:self.mid_idx,i], fmt='r.')
                plt.errorbar(self.V[self.mid_idx:], n[self.mid_idx:, i] + i, yerr=nErr[self.mid_idx:,i],fmt='b.')
                plt.xlabel('Voltage')
                plt.ylabel('Occupation')
            plt.figure(figsize=FIGSIZE)
            for i in range(len(Q[0])):
                plt.errorbar(self.V[:self.mid_idx], -Q[:self.mid_idx, i] + i, yerr=QErr[:self.mid_idx, i], fmt='r.')
                plt.errorbar(self.V[self.mid_idx:], -Q[self.mid_idx:, i] + i, yerr=QErr[self.mid_idx:, i], fmt='b.')
                plt.xlabel('Voltage')
                plt.ylabel('Chagre')


            if self.full_I is not None:
                plt.figure(figsize=FIGSIZE)
                for i in range(len(self.full_I)):
                    plt.plot(self.V[:self.mid_idx], self.full_I[i,:self.mid_idx], 'o',
                            self.V[self.mid_idx:], self.full_I[i,self.mid_idx:], '*')
                    plt.xlabel('Voltage')
                    plt.ylabel('Current')

    def power_law_fit(self, I=None, IErr=None, V=None, mid_idx=None, err=True, errorevery=10, Vlabel='Voltage', Ilabel='Current', shift=0,
                        fmt_up='r.', fmt_down='b.', plot=False):
        if I is None:
            I = self.I
        if IErr is None:
            IErr = self.IErr
        if V is None:
            V = self.V
        if mid_idx is None:
            mid_idx = self.mid_idx
        VthresholdDown = self.calc_threshold_voltage_down(I, IErr, V, mid_idx)[0]
        VthresholdUp = self.calc_threshold_voltage_up(I, IErr, V, mid_idx)[0]
        if VthresholdUp == 0 or VthresholdDown == 0:
            return [0,0],[0,0],[0,0],[0,0]
        VUp = V[:mid_idx]
        IUp = I[:mid_idx]
        IErrUp = IErr[:mid_idx]

        IUp = IUp[VUp > VthresholdUp]
        IErrUp = IErrUp[VUp > VthresholdUp]
        VUp = (VUp[VUp > VthresholdUp] - VthresholdUp) / VthresholdUp
        VDown = V[mid_idx:]
        IDown = I[mid_idx:]
        IErrDown = IErr[mid_idx:]

        IDown = IDown[VDown > VthresholdDown]
        IErrDown = IErrDown[VDown > VthresholdDown]
        VDown = (VDown[VDown > VthresholdDown] - VthresholdDown) / VthresholdDown
        if VDown.size < 30 or VUp.size < 30:
            return [0,0],[0,0],[0,0],[0,0]
        VUp_for_fit = VUp[:50]
        VDown_for_fit = VDown[-50:]
        IUp_for_fit = IUp[:50] - IUp[0]
        IDown_for_fit = IDown[-50:] - IDown[-1]

        def power_law(x, a, b):
            return a * x ** b
        try:
            paramsUp, covUp = curve_fit(f=power_law, xdata=VUp_for_fit, ydata=IUp_for_fit,
                                        p0=[0, 0], bounds=(-np.inf, np.inf))
            paramsDown, covDown = curve_fit(f=power_law, xdata=VDown_for_fit, ydata=IDown_for_fit,
                                            p0=[0, 0], bounds=(-np.inf, np.inf))
        except RuntimeError:
            return [0,0],[0,0],[0,0],[0,0]
        if plot:
            if err:
                plt.errorbar(VUp, IUp + shift, fmt=fmt_up, yerr=IErrUp,
                             errorevery=errorevery, label="up")
                plt.errorbar(VDown, IDown + shift, fmt=fmt_down,
                             yerr=IErrDown,
                             errorevery=errorevery, label="down")
            else:
                plt.plot(VUp, IUp + shift, fmt_up, label="up")
                plt.plot(VDown, IDown + shift, fmt_down, label="down")
            print("Up fit: {}x^{}".format(*paramsUp))
            plt.plot(VUp_for_fit, power_law(VUp_for_fit, *paramsUp) + IUp[0], label="fit up")
            print("Down fit: {}x^{}".format(*paramsDown))
            plt.plot(VDown_for_fit, power_law(VDown_for_fit, *paramsDown) + IDown[-1], label="fit down")
            plt.xlabel(Vlabel)
            plt.ylabel(Ilabel)
        return paramsUp, np.sqrt(np.diag(covUp)), paramsDown, np.sqrt(np.diag(covDown))

    def calc_param_average(self, name):
        if name in ["C","R"]:
            combined = np.hstack((np.array(self.get_array_param(name + "h")).flatten(),np.array(self.get_array_param(name + "v")).flatten()))
            avg = np.average(combined)
        else:
            avg = np.average(self.get_array_param(name))
        return avg

    def calc_paths(self):
        if self.Imaps is None:
            print("Imaps wasn't loaded")
            return False
        path_dict = dict()
        if len(self.Imaps.shape) < 3:
            self.Imaps = self.Imaps.reshape((1,1,self.Imaps.size))
        for v,Imap in zip(self.V, self.Imaps):
            paths, paths_current = findPaths(Imap, self.rows, self.columns)
            for path, path_current in zip(paths, paths_current):
                path_tup = tuple(path)
                current = np.min(path_current)
                if path_tup not in path_dict:
                    path_dict[path_tup] = [[],[]]
                path_dict[path_tup][0].append(v)
                path_dict[path_tup][1].append(current)
        return self.filter_paths(path_dict, 30)

    def filter_paths(self, path_dict, percentage_threshold):
        currents = []
        if len(path_dict) > 0:
            for path in path_dict:
                currents.append(np.max(path_dict[path][1]))
            eps = percentage_threshold * np.max(currents) / 100
            to_delete = []
            for path in path_dict:
                if np.max(path_dict[path][1]) < eps:
                    to_delete.append(path)
            for path in to_delete:
                del path_dict[path]
        return path_dict

    def plot_current_by_paths(self, fig=None, ax=None):
        if self.path_dict is None:
            self.path_dict = self.calc_paths()
        if fig is None:
            fig = plt.figure(figsize=FIGSIZE)
        if ax is None:
            ax = fig.add_subplot(111)
        colors = matplotlib.cm.gist_rainbow(np.linspace(0, 1, len(self.path_dict.keys())))
        for idx, path in enumerate(sorted(self.path_dict.keys())):
            v = np.array(self.path_dict[path][0])
            i = np.array(self.path_dict[path][1])
            mid_idx = np.argmax(v)
            ax.plot(v[:mid_idx], i[:mid_idx], linestyle='dotted', color=colors[idx])
            ax.plot(v[mid_idx:], i[mid_idx:], linestyle='dotted', color=colors[idx])
        return fig, ax

class MultiResultAnalyzer:
    """ Used for statistical analysis of results from many simulations"""
    def __init__(self, directories_list, files_list, relevant_running_params=None, relevant_array_params=None, out_directory = None,
                 groups=None, group_names=None, resistance_line=0, full=False, graph=False, reAnalyze=False, filter=lambda x:True):
        """
        Initializing analyzer
        :param directories_list: list of directories for result files
        :param files_list: list of result names
        :param relevant_array_params: list of array parameters that are relevant
        :param relevant_running_params: list of array parameters that are relevant
        :param out_directory: directory for results
        """
        self.outDir = out_directory
        if out_directory is not None:
            if not os.path.isdir(out_directory):
                os.mkdir(out_directory)
        self.directories = directories_list
        self.fileNames = files_list
        self.scores = dict()
        for score in SCORE_NAMES:
            self.scores[score] = {SCORE_VAL:[],SCORE_HIGH_ERR:[],SCORE_LOW_ERROR:[]}
        if relevant_array_params:
            self.arrayParams = {p: [] for p in relevant_array_params}
        else:
            self.arrayParams = dict()
        if relevant_running_params:
            self.runningParams = {p: [] for p in relevant_running_params}
        else:
            self.runningParams = dict()
        self.disorders = {p: [] for p in ["C","R","CG","VG"]}
        self.averages = {p: [] for p in ["C", "R", "CG", "VG"]}
        self.groups = np.array(groups)
        self.groupNames = group_names
        self.small_v_fit_up = []
        self.small_v_cov_up = []
        self.small_v_fit_down = []
        self.small_v_cov_down = []
        load_params = (relevant_array_params is not None) or (relevant_running_params is not None)
        self.load_data(resistance_line, full=full, graph=graph, load_params=load_params, reAnalyze=reAnalyze,
                       filter=filter)

    def load_data(self, line=0, full=False, graph=False, load_params=True, reAnalyze=False, filter=lambda x:True):
        """
        Loading results and calculating scores
        """
        missing_idx = []
        idx = 0
        for directory, fileName in zip(self.directories, self.fileNames):
            try:
                processor = SingleResultsProcessor(directory, fileName, fullOutput=full, graph=graph, reAnalyze=reAnalyze)
                processor.load_results(reAnalyze=reAnalyze)
                if load_params:
                    processor.load_params()
            except MissingFilesException:
                missing_idx.append(idx)
                continue
            if not filter(processor):
                missing_idx.append(idx)
                continue
            for score in self.scores:
                val, high, low = processor.calc_score(score)
                self.scores[score][SCORE_VAL].append(val)
                self.scores[score][SCORE_HIGH_ERR].append(high)
                self.scores[score][SCORE_LOW_ERROR].append(low)


            if load_params:
                for p in self.arrayParams:
                    self.arrayParams[p].append(processor.get_array_param(p))
                for p in self.runningParams:
                    self.runningParams[p].append(processor.get_running_param(p))
                for p in self.disorders:
                    std, avg = self.calc_disorder_and_average(p, processor)
                    self.disorders[p].append(std)
                    self.averages[p].append(avg)
            idx += 1
        for i in missing_idx:
            del self.directories[i]
            del self.fileNames[i]
    def calc_disorder_and_average(self, name, processor):
        if name in ["C","R"]:
            combined = np.hstack((np.array(processor.get_array_param(name + "h")).flatten(),np.array(processor.get_array_param(name + "v")).flatten()))
            std = np.std(combined)
            avg = np.average(combined)
        else:
            std = np.std(processor.get_array_param(name))
            avg = np.average(processor.get_array_param(name))
        return std, avg

    def get_relevant_scores(self, scoreName):
        if scoreName in self.scores:
            score = self.scores[scoreName]
            return score[SCORE_VAL], score[SCORE_HIGH_ERR], score[SCORE_LOW_ERROR]
        else:
            raise NotImplementedError



    def plot_hystogram(self, scoreName, runinng_parameter_vals, array_parameter_vals, title):
        """
        Plotting hystogram of the given score for results with parameters matching the given values
        :param scoreName: name of score (hysteresis, jump, blockade)
        :param (running/array)_parameter_vals: dictionary with {parameter_name:list_of_parameter_values}.
        For non-given parameters all values will be considered as valid
        """
        score= self.get_relevant_scores(scoreName)[SCORE_VAL]
        indices = set(range(score.size))
        for p in runinng_parameter_vals:
            good_inds = set()
            good_vals = runinng_parameter_vals[p]
            for ind,val in enumerate(self.runningParams[p]):
                if val in good_vals:
                    good_inds.add(ind)
            indices = indices.intersection(good_inds)
        for p in array_parameter_vals:
            good_inds = set()
            good_vals = array_parameter_vals[p]
            for ind,val in enumerate(self.arrayParams[p]):
                if val in good_vals:
                    good_inds.add(ind)
            indices = indices.intersection(good_inds)
        indices = list(indices)
        scores = score[indices]
        bins = max(len(scores)//10, 5)
        fig = plt.figure()
        plt.hist(scores, bins=bins)
        plt.title(title)
        plt.xlabel("score")
        plt.ylabel("freq.")
        plt.savefig(os.path.join(self.outDir, title.replace(' ', '_') + '.png'))
        plt.close(fig)
        np.save(os.path.join(self.outDir,title), np.array(scores))

    def plot_score_by_groups(self, score1_name, score2_name):
        score1 = self.get_relevant_scores(score1_name)[SCORE_VAL]
        score2 = self.get_relevant_scores(score2_name)[SCORE_VAL]
        plt.figure()
        for i in range(np.max(self.groups)+1):
            plt.scatter(score1[self.groups == i],score2[self.groups == i], label=self.groupNames[i])
        plt.legend()
        plt.xlabel(score1_name)
        plt.ylabel(score2_name)
        plt.savefig(os.path.join(self.outDir, score1_name + "_vs_" + score2_name + '.png'))



    def get_y_label(self, score_name):
        if score_name == "hysteresisArea":
            return "hysteresis loop area"
        elif score_name == "jump ":
            return "voltage diff for biggest jump"
        elif score_name == "blockade":
            return "threshold voltage"


    def plot_score(self, score_name, xs=None, xlabel=None, ylabel=None, label="", shift=0, fmt='.'):
        y, y_high_err, y_low_err = self.get_relevant_scores(score_name)
        yerror = np.zeros((2,len(y)))
        yerror[0,:] = y_low_err
        yerror[1,:] = y_high_err
        x = xs if xs is not None else np.arange(len(y))
        y = np.array(y)
        x = np.array(x)
        yerror = yerror[:,y>0]
        x = x[y>0]
        y = y[y > 0]
        if ylabel is None:
            ylabel = score_name
        if xlabel is None:
            xlabel = "example num."
        plt.errorbar(x, np.array(y) + shift, yerr=yerror, fmt=fmt, label=label)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)

    def plot_scores_together(self, score1, score2, average_results=True,  fig=None, ax=None, window_size=0.01,
                             fmt='.'):
        if ax is None:
            fig,ax = plt.subplots(figsize=FIGSIZE)
        x, x_high_err, x_low_err = self.get_relevant_scores(score1)
        y, y_high_err, y_low_err = self.get_relevant_scores(score2)
        yerror = np.zeros((2, len(y)))
        yerror[0, :] = y_low_err
        yerror[1, :] = y_high_err
        xerror = np.zeros((2, len(x)))
        xerror[0, :] = x_low_err
        xerror[1, :] = x_high_err
        x = np.array(x)
        y = np.array(y)
        relevant = np.logical_and(x>0, y>0)
        xerror = xerror[:, relevant]
        yerror = yerror[:, relevant]
        x = x[relevant]
        y = y[relevant]
        if average_results:
            x, y, yerror = self.running_average(x, y, yerror, window_size=window_size)
        ax.errorbar(x, y, yerr=yerror, fmt=fmt)
        return ax, fig

    def plot_score_by_parameter(self, score_name, parameter_names, runningParam=True,
                                average_results=True, normalizing_parameter_names=None, fig=None, ax=None):
        """
        Plots score as a function of given parameter
        :param score_name: relevant score name (hysteresis, jump, blockade)
        :param parameter_name: relevant parameter
        """
        self.print_max_score_file(score_name)
        if len(parameter_names) == 1:
            if ax is None:
                fig, ax = plt.subplots(figsize=FIGSIZE)
            parameter_name = parameter_names[0]
            y, y_high_err, y_low_err = self.get_relevant_scores(score_name)
            if runningParam:
                normalization = np.array(self.runningParams[normalizing_parameter_names[0]]) if\
                    normalizing_parameter_names is not None else 1
                x = np.array(self.runningParams[parameter_name])/normalization
            else:
                normalization = np.array(self.arrayParams[normalizing_parameter_names[0]]) if \
                    normalizing_parameter_names is not None else 1
                x = np.array(self.arrayParams[parameter_name])/normalization
            if average_results:
                x, y, y_high_err, y_low_err = self.average_similar_results(x, y, y_high_err, y_low_err)
                x = np.array(x)
            y = np.array(y)
            yerror = np.zeros((2, len(y)))
            yerror[0, :] = y_low_err
            yerror[1, :] = y_high_err
            yerror = yerror[:, y > 0]
            x = x[y > 0]
            y = y[y > 0]
            plt.errorbar(x, y, yerr=yerror, fmt='.',)
        elif len(parameter_names) == 2:
            if ax is None:
                fig = plt.figure(figsize=FIGSIZE)
                ax = fig.add_subplot(111, projection="3d")
            z, z_high_err, z_low_err = self.get_relevant_scores(score_name)
            if runningParam:
                x = np.array(self.runningParams[parameter_names[0]])
                y = np.array(self.runningParams[parameter_names[1]])
            else:
                x = np.array(self.arrayParams[parameter_names[0]])
                y = np.array(self.arrayParams[parameter_names[1]])
            points = [(x[i],y[i]) for i in range(len(x))]
            if average_results:
                points, z, z_high_err, z_low_err = self.average_similar_results(points, z, z_high_err, z_low_err)
            zerror = np.zeros((2, len(z)))
            zerror[0, :] = z_low_err
            zerror[1, :] = z_high_err
            z = np.array(z)
            x = np.array([point[0] for point in points])
            y = np.array([point[1] for point in points])
            x = x[z>0]
            y = y[z>0]
            zerror = zerror[:,z>0]
            z = z[z>0]
            plt.xlabel(parameter_names[0])
            plt.ylabel(parameter_names[1])

        else:
            raise NotImplementedError
        return ax, fig

    def plot_results_by_disorder(self, parameters, scores, plot3D=False, fig=None, ax=None, average=True,
                                 window_size=None,fmt='.'):
        for score in scores:
            self.print_max_score_file(score)
        norm = []
        for param in parameters:
            if param in ("R", "C"):
                norm.append(np.array(self.averages[param]))
            elif param == "CG":
                norm.append(np.array(self.averages["C"]))
            else:
                norm.append(1)
        if fig is None:
            fig = plt.figure(figsize=FIGSIZE)
        if plot3D and len(parameters) ==2:
            if ax is None:
                ax = fig.add_subplot(111, projection='3d')
            for score in scores:
                x = np.array(self.disorders[parameters[0]]) / norm[0]
                y = np.array(self.disorders[parameters[1]]) / norm[1]
                z, z_high_err, z_low_err = self.get_relevant_scores(score)
                z = np.array(z)
                zerror = np.zeros((2, len(z)))
                zerror[0, :] = z_low_err
                zerror[1, :] = z_high_err
                x = x[z>0]
                y = y[z>0]
                zerror = zerror[:,z>0]
                z = z[z>0]
                if average:
                    if window_size is None:
                        window_size = (0,0)
                    x,y,z,zerror = self.running_average_3d(x,y,z,zerror,window_size)
                ax.scatter(x, y, z)
                ax.set_ylim(0,1.2)
                ax.set_xlim(0,0.4)
        else:
            if ax is None:
                ax = fig.add_subplot(111)
            for param, normalization in zip(parameters,norm):
                for score in scores:
                    x = np.array(self.disorders[param])/normalization
                    y, y_high_err, y_low_err = self.get_relevant_scores(score)
                    y = np.array(y)
                    yerror = np.zeros((2, len(y)))
                    yerror[0, :] = y_low_err
                    yerror[1, :] = y_high_err
                    x = x[y>0]
                    yerror = yerror[:,y>0]
                    y = y[y>0]
                    if average:
                        if window_size is None:
                            window_size = 0
                        x, y, yerror = self.running_average(x, y, yerror, window_size)

                    ax.errorbar(x,y, yerr=yerror, fmt=fmt, markersize=14)
        return ax, fig

    def plot_results_by_average(self, parameters, scores, plot3D=False, fig=None, ax=None, average=True,
                                window_size=None, fmt='.'):
        for score in scores:
            self.print_max_score_file(score)
        norm = []
        for param in parameters:
            if param == "CG":
                norm.append(np.array(self.averages["C"]))
            else:
                norm.append(1)
        if fig is None:
            fig = plt.figure(figsize=FIGSIZE)
        if plot3D and len(parameters) ==2:
            if ax is None:
                ax = fig.add_subplot(111, projection='3d')
            for score in scores:
                x = np.array(self.averages[parameters[0]])/norm[0]
                y = np.array(self.averages[parameters[1]])/norm[1]
                z, z_high_err, z_low_err = self.get_relevant_scores(score)
                z = np.array(z)
                zerror = np.zeros((2, len(z)))
                zerror[0, :] = z_low_err
                zerror[1, :] = z_high_err
                x = x[z>0]
                y = y[z>0]
                zerror = zerror[:, z > 0]
                z = z[z>0]
                if average:
                    if window_size is None:
                        window_size = (0,0)
                    x,y,z,zerror = self.running_average_3d(x,y,z,zerror,window_size)
                ax.scatter(x, y, z)
        else:
            if ax is None:
                ax = fig.add_subplot(111)
            for param,normalization in zip(parameters, norm):
                for score in scores:
                    x = np.array(self.averages[param])/normalization
                    y, y_high_err, y_low_err = self.get_relevant_scores(score)
                    y = np.array(y)
                    yerror = np.zeros((2, len(y)))
                    yerror[0, :] = y_low_err
                    yerror[1, :] = y_high_err
                    x = x[y > 0]
                    yerror = yerror[:, y > 0]
                    y = y[y > 0]
                    if average:
                        if window_size is None:
                            window_size = 0
                        x, y, yerror = self.running_average(x,y,yerror,window_size)
                    ax.errorbar(x,y, yerr=yerror, fmt=fmt)
        return ax, fig

    def print_max_score_file(self, score):
        index = np.argmax(self.scores[score]["score"])
        val = np.max(self.scores[score]["score"])
        print("File with maximum " + score + " is:")
        print(os.path.join(self.directories[index],self.fileNames[index]))

    def average_similar_results(self, xs, ys, ys_high_err, ys_low_err):
        ys = np.array(ys)
        ys_high_err = np.array(ys_high_err)
        ys_low_err = np.array(ys_low_err)
        new_x = []
        new_y = []
        new_y_high_err = []
        new_y_low_err = []
        for x in xs:
            if x not in new_x:
                new_x.append(x)
                indices = [i for i, val in enumerate(xs) if val == x and ys[i] != 0]
                relevant_ys = ys[indices]
                relevant_high_err = ys_high_err[indices]
                relevant_low_err = ys_low_err[indices]
                new_y.append(np.average(relevant_ys))
                new_y_high_err.append(np.sqrt(np.average(relevant_high_err**2)/relevant_high_err.size))
                new_y_low_err.append(np.sqrt(np.average(relevant_low_err**2)/relevant_low_err.size))
        return new_x, new_y, new_y_high_err, new_y_low_err

    def running_average(self, x, y, yErr, window_size=0):
        if x.size ==0:
            return x, y, yErr
        x_min = np.min(x)
        x_max = np.max(x)
        if window_size == 0:
            window_size = (x_max-x_min)/10
        low = x_min
        high = x_min + window_size
        new_xs = []
        new_ys = []
        new_ysErr = []
        while high <= x_max:
            relevant = np.logical_and(x >= low, x <= high)
            if x[relevant].size > 3:
                new_xs.append(np.average(x[relevant]))
                new_ys.append(np.average(y[relevant]))
                new_ysErr.append(np.std(y[relevant]/len(y[relevant])))
            low += window_size
            high += window_size
        return np.array(new_xs), np.array(new_ys), np.array(new_ysErr)

    def running_average_3d(self, x, y, z, zErr, window_size=(0,0)):
        zErr = np.max(zErr, axis=0)
        x_min = np.min(x)
        x_max = np.max(x)
        y_min = np.min(y)
        y_max = np.max(y)
        if window_size == (0,0):
            window_size = ((x_max-x_min)/10, (y_max-y_min)/10)
        low_x = x_min
        high_x = x_min + window_size[0]
        low_y = y_min
        high_y = y_min + window_size[1]
        new_xs = []
        new_ys = []
        new_zs = []
        new_zsErr = []
        while high_x <= x_max:
            while high_y <= y_max:
                relevant = np.logical_and(np.logical_and(y>=low_y, y <= high_y),np.logical_and(x>=low_x,x<=high_x))
                new_xs.append(np.average(x[relevant]))
                new_ys.append(np.average(y[relevant]))
                new_z_variance = 1/np.sum(1/zErr[relevant]**2)
                new_zs.append(new_z_variance*np.sum(z[relevant]/zErr[relevant]**2))
                new_zsErr.append(np.sqrt(new_z_variance))
        return np.array(new_xs), np.array(new_ys), np.array(new_zs), np.array(new_zsErr)


def get_filter(param_names, param_vals, score_names, minimum_score_vals):
    def filter(processor):
        for param_name, param_val in zip(param_names, param_vals):
            if processor.get_running_param(param_name) != param_val:
                return False
        for score_name, score_val in zip(score_names, minimum_score_vals):
            score, _, _ = processor.calc_score(score_name)
            if score < score_val:
                return False
        return True
    return filter


def resistance_temperature_analysis(directory, file_names):
    I = []
    T = []
    V = []
    Vthreshold = []
    for file in file_names:
        s = SingleResultsProcessor(directory, file)
        I.append(s.I)
        T.append(s.get_running_param("T"))
        V.append(s.V)
        temp = s.V[s.I>0]
        Vthreshold.append(temp[0])
    I = np.array(I)
    V = np.array(V)
    T = np.array(T)
    for i in range(0,I.shape[1]//4,50):
        It = I[:,i]
        Vt = V[:,i]
        plt.semilogy(1/T, np.abs(Vt/It))

def findPaths(Imap, rows, columns, eps=0.001):
    paths = []
    paths_current = []
    def array_DFS(i ,j , path, path_current):
        if j == columns:
            path_copy = path.copy()
            path_current_copy = path_current.copy()
            path_copy.append((i,j))
            paths.append(path_copy)
            paths_current.append(path_current_copy)
            path_flow = np.min(path_current_copy)
            first_row = path_copy[0][0]
            Imap[first_row*2,0] -= path_flow
            for node_idx in range(len(path_copy)-1):
                flow_row, flow_column = path_copy[node_idx]
                flow_row += path_copy[node_idx+1][0]
                flow_column += max(path_copy[node_idx+1][1]-flow_column,0)
                if Imap[flow_row, flow_column] > 0:
                    Imap[flow_row, flow_column] -= path_flow
                else:
                    Imap[flow_row, flow_column] += path_flow
            return True
        elif (i,j) in path or i < 0 or j < 0 or i == rows:
            return False
        else:
            path.append((i,j))
            if Imap[i*2,j] < -eps:
                path_current.append(-Imap[i*2 , j ])
                if array_DFS(i,j-1,path, path_current):
                    return True
                del path_current[-1]
            if i > 0 and Imap[i*2-1,j] < -eps:
                path_current.append(-Imap[i*2 - 1, j ])
                if array_DFS(i - 1, j, path, path_current):
                    return True
                del path_current[-1]
            if i < rows - 1 and Imap[i*2 + 1,j] > eps:
                path_current.append(Imap[i*2 + 1, j])
                if array_DFS(i+1, j, path, path_current):
                    return True
                del path_current[-1]
            if Imap[i * 2, j + 1] > eps:
                path_current.append(Imap[i * 2, j + 1])
                if array_DFS(i, j + 1, path, path_current):
                    return True
                del path_current[-1]
        del path[-1]
        return False
    for i in range(rows):
        if Imap[i*2, 0] > eps:
            path = []
            path_current = [Imap[i*2, 0]]
            while array_DFS(i, 0, path, path_current) and Imap[i*2, 0] > eps:
                path = []
                path_current = [Imap[i * 2, 0]]
    return paths, paths_current

def plot_paths(paths, ax, colors):
    shift = 0
    for path_idx,path in enumerate(paths):
        i_start = path[0][0]
        ax.arrow(6*shift, 6 * (i_start + shift), 4, 0,
                 head_width=1, head_length=0.5, color=colors[path_idx], width=0.1)
        for idx in range(len(path)-1):
            i1 = path[idx][0]
            j1 = path[idx][1]
            i2 = path[idx+1][0]
            j2 = path[idx+1][1]
            ax.arrow(4+6*(j1+ shift), 6*(i1+shift), 6*(j2-j1), 6*(i2-i1),
                     head_width=1, head_length=0.5, color=colors[path_idx], width=0.1)
        shift += 0.02
    return ax

def average_diff(x, window_size=5):
    if window_size < 5:
        window_size = 5
    if 2*window_size > len(x):
        window_size = len(x)//2
    if window_size == 0:
        return np.diff(x)

    def preper_kernel(x, window_size):
        kernel = (1 / window_size) * np.ones((2 * window_size,))
        kernel[window_size:] *= -1
        padded = np.pad(x, window_size - 1, 'edge')
        return padded, kernel

    padded, kernel = preper_kernel(x, window_size)
    result = np.convolve(padded, kernel, mode='valid')
    return result




def bistabilityAnalysis(x):
    avg = np.average(x)
    var = np.var(x)
    dist_from_mean = (x - avg)**2
    bins = [x]
    if np.max(dist_from_mean) > var:
        single_pick_model = GaussianMixture(n_components=1)
        double_pick_model = GaussianMixture(n_components=2)
        single_pick_model.fit_predict(x.reshape(-1, 1))
        double_pick_labels = double_pick_model.fit_predict(x.reshape(-1, 1))
        if double_pick_model.aic(x.reshape(-1, 1)) < single_pick_model.aic(x.reshape(-1, 1)):
            bins = [x[double_pick_labels==0], x[double_pick_labels==1]]
            if len(bins[0]) > 1 and len(bins[1]) > 1:
                max_std = max(np.std(bins[0]), np.std(bins[1]))
                mean_dist = np.abs(np.average(bins[0]) - np.average(bins[1]))
                if max_std > mean_dist:
                    bins = [x]
        else:
            bins = [x]
    return bins


def approx_jumps_separation(R1, R2, C1, C2, CG):
    if R2>R1:
        return 1 / ((CG / 2) * ((R2 - R1) / (R2 + R1)) + C2)
    else:
        return 1 / ((CG / 2) * ((R1 - R2) / (R2 + R1)) + C1)


def approx_threshold_voltage(R1, R2, C1, C2, CG, VG, nini):
    q1 = -CG * VG + 0.5 * (C1 + C2 + CG) / (C1 + C2) - nini
    q2 = CG * VG + 0.5 * (C1 + C2 + CG) / (C1 + C2) + nini
    Vt1 = q1 / (C1 + CG / 2)
    Vt2 = q2 / (C2 + CG / 2)
    if Vt2 < Vt1:
        if R1 > R2:
            return Vt2
        else:
            nini = 1
            q1 = -CG * VG + 0.5 * (C1 + C2 + CG) / (C1 + C2) - nini
            q2 = CG * VG + 0.5 * (C1 + C2 + CG) / (C1 + C2) + nini
            Vt1 = q1 / (C1 + CG / 2)
            Vt2 = q2 / (C2 + CG / 2)
            return min(Vt1, Vt2)
    else:
        if R2 > R1:
            return Vt1
        else:
            nini = -1
            q1 = -CG * VG + 0.5 * (C1 + C2 + CG) / (C1 + C2) - nini
            q2 = CG * VG + 0.5 * (C1 + C2 + CG) / (C1 + C2) + nini
            Vt1 = q1 / (C1 + CG / 2)
            Vt2 = q2 / (C2 + CG / 2)
            return min(Vt1, Vt2)

def approx_jumps_height(R1, R2, C1, C2, CG):
    return 1 / ((C1 + C2 + CG) * (R1 + R2))

def approx_jump_up_down(R1, R2, C1, C2, CG, V):
    Rmax = max(R1, R2)
    Rmin = min(R1, R2)
    Ugap = np.abs((Rmax/(C1+C2)+ Rmin * V)/(Rmax-Rmin) - (C1+C2+2*CG)*(np.sqrt((R1*R2*V)/(CG*(C1+C2+CG)*(C1+C2)*(R1-R2)**2))))
    # Ugap = ((Rmax/(C1+C2)+ ((C1+C2+CG)/CG)*Rmin * V - 2*(np.sqrt((R1*R2*V*(C1+C2+CG))/(CG*(C1+C2))))))/(Rmax-Rmin)
    return CG * Ugap * approx_jumps_separation(R1, R2, C1, C2, CG)

def approx_hysteresis_total_area(R1, R2, C1, C2, CG, VG, Vmax):
    Vmax = min((max(R1, R2)*CG)/(min(R1, R2)*(C1+C2)*(C1+C2+CG)),Vmax)
    Vmin = 1.1
    dV = approx_jumps_separation(R1,R2,C1,C2,CG)
    V = Vmin
    area = 0
    while V < Vmax:
        area += approx_jump_up_down(R1,R2,C1,C2,CG,V)*approx_jumps_height(R1,R2,C1,C2,CG)
        V += dV
    return area

def significant_figures(x, xerr):
    factor = 1
    if xerr > 5:
        while xerr*factor > 5:
            factor /= 10
    elif xerr < 0.5:
        while xerr*factor < 0.5:
            factor *=10
    return int(np.round(x*factor))/factor, int(np.round(xerr*factor))/factor


def getOptions():
    parser = OptionParser(usage= "usage: %prog [options]")
    parser.add_option("-o", "--output-folder", dest="output_folder",
                      help="Output folder [default: no outputcurrent folder]",
                      default='')
    parser.add_option("-i", "--input", dest="directories", action="append",
                      help="input directory, can give more than one",
                      default=[])
    parser.add_option("--file", dest="file_names", action="append", default=[],
                      help="List of input file names [default: all files in input directory]")
    parser.add_option("--action", dest="action", help="The processing action, one of: " + str(ACTIONS))
    parser.add_option("--score", dest="scores",
                      help="The score to plot, can give more than one. For score plotting actions. One of: " + str(SCORE_NAMES),
                      action="append", default=[])
    parser.add_option("--parameter", dest="parameters", action="append", default=[],
                      help="Parameter  or parameters for x axis of plot. For score plotting actions")
    parser.add_option("--normalizing", dest="normalizing", action="append", default=[],
                      help="normalize the plot parameters by these parameteres, given in the same order as parameters.")
    parser.add_option("--plot3d", dest="plot3d", action="store_true", default=False,
                      help="Plots score in 3d plot by 2 parameters. For score plotting actions")
    parser.add_option("--full", dest="full", action="store_true", default=False,
                      help="If true full output will be loaded (including occupation and gate charge of each island")
    parser.add_option("--average", dest="average", action="store_true", default=False,
                      help="If true similar results would be averaged")
    parser.add_option("--by-occupation", dest="by_occupation", action="store_true", default=False,
                      help="Only used for plotting jumps. If ture jumps would be plotted by island occupation,"
                           " otherwise by paths current.")
    parser.add_option("--by-path-occupation", dest="by_path_occupation", action="store_true", default=False,
                      help="Only used for plotting jumps. If ture jumps would be plotted by path occupation,"
                           " otherwise by paths current.")
    parser.add_option("--re-analyze", dest="re_analyze", action="store_true", default=False,
                      help="If true performs multigaussian reanalysis result for the IV curve")
    parser.add_option("--save-re-analysis", dest="save_re_analysis", action="store_true", default=False,
                      help="If true saves the multigaussian reanalysis result for the IV curve")
    parser.add_option("--window-size", dest="window", type=float, default=0,
                      help="Window size for running average implementation. Relevant for plotting score by disorder "
                           "or average actions. [default: 1/10 from the range of the given data")
    parser.add_option("--filter-param", dest="filter_param", default=[], action="append",
                      help="Parameter by which to filter the results (only results where parameter equals the given "
                           "value would be processed")
    parser.add_option("--filter-val", dest="filter_val", default=[], action="append", type=float,
                      help="Parameter value by which to filter the results (only results where parameter equals the"
                           " given value would be processed")
    parser.add_option("--filter-score", dest="filter_score", default=[], action="append",
                      help="Score by which to filter the results (only results where score is bigger than the given "
                           "value would be processed")
    parser.add_option("--min-score-val", dest="min_score_val", default=[], action="append", type=float,
                      help="Score value by which to filter the results (only results where score is bigger than the"
                           " given value would be processed")
    parser.add_option("--xlabel", dest="xlabel", default="",
                      help="xlabel for plot")
    parser.add_option("--ylabel", dest="ylabel", default="",
                      help="ylabel for plot")
    parser.add_option("--zlabel", dest="zlabel", default="",
                      help="zlabel for plot")
    return parser.parse_args()




if __name__ == "__main__":
    options, args = getOptions()
    if options.output_folder:
        if os.path.exists(options.output_folder):
           if not os.path.isdir(options.output_folder):
               print("Can't create output directory, a file with the same name already exists.")
               exit(0)
        else:
            os.mkdir(options.output_folder)
    directories = options.directories
    inner_directories = [os.path.join(directory,name) for directory in directories for name in os.listdir(directory) if os.path.isdir(os.path.join(directory,name))]
    directories += inner_directories
    action = options.action
    if len(options.file_names) > 0:
        file_names = [options.file_names]
    else:
        if action == "plot_RT":
            file_names = [[name.replace("_IT.png", "") for name in os.listdir(directory) if "_IT.png" in name] for
                          directory in directories]
        else:
            file_names = [[name.replace("_IV.png","") for name in os.listdir(directory) if "_IV.png" in name] for directory in directories]
    filter = get_filter(options.filter_param, options.filter_val, options.filter_score, options.min_score_val)
    if action == "perpendicular_current":
        full = True
        for directory, names in zip(directories, file_names):
            for name in names:
                s = SingleResultsProcessor(directory, name, fullOutput=full, vertCurrent=True)
                if filter(s):
                    s.plot_resistance(out=os.path.join(directory,name))
    elif action == "plot_RT":
        directory = "/home/kasirershahar/University/Research/simulation_results/single_island/single_island_IT"
        files = [["array_1_1_cg_1_v_0.1","array_1_1_cg_1_v_0.1_r_disorder","array_1_1_cg_10_v_0.1",
                  "array_1_1_cg_10_v_0.1_r_disorder"],
                 ["array_1_1_cg_1_v_0.5", "array_1_1_cg_1_v_0.5_r_disorder", "array_1_1_cg_10_v_0.5",
                  "array_1_1_cg_10_v_0.5_r_disorder"],
                 ["array_1_1_cg_1_v_1", "array_1_1_cg_1_v_1_r_disorder", "array_1_1_cg_10_v_1",
                  "array_1_1_cg_10_v_1_r_disorder"]]
        fig, ax = plt.subplots(3,1, figsize=FIGSIZE)
        shift = 1
        shift_factor = [700,700,100]
        labels =["$\\frac{C_G}{C_1+C_2} = 1$ $\\frac{R_2}{R_1}=1$", "$\\frac{C_G}{C_1+C_2} = 1$ $\\frac{R_2}{R_1}=9$", "$\\frac{C_G}{C_1+C_2} = 10$ $\\frac{R_2}{R_1}=1$",
                 "$\\frac{C_G}{C_1+C_2} = 10$ $\\frac{R_2}{R_1}=9$"]
        colors = matplotlib.cm.gist_rainbow(np.linspace(0, 1, 4))
        for index in range(3):
            for i,name in enumerate(files[index]):
                s = SingleResultsProcessor(directory, name, fullOutput=options.full, vertCurrent=False, graph=False,
                                           reAnalyze=options.re_analyze, IT=True)
                if filter(s):
                    v = s.plot_RT(ax[index], shift=shift, err=True, color=colors[i])
                    shift *= shift_factor[index]
                ax[index].set_title("$V\\frac{C_1+C_2}{e} = %1.1f$" %v, position=(0.5, 0.75))
                ax[index].set_xlabel("$T\\frac{C_1+C_2}{e^2}$")
                ax[index].set_ylabel("$\\frac{R}{R_1+R_2}$")
            index += 1
            shift=1
        ax[0].set_xlabel('')
        ax[0].set_xticks([])
        ax[0].set_xlim(0,0.3)
        ax[0].set_ylim(1,1e15)
        ax[1].set_xlabel('')
        ax[1].set_xticks([])
        ax[1].set_xlim(-0.003, 0.3)
        ax[1].set_ylim(1,1e15)
        ax[2].set_ylim(1,1e10)
        ax[2].set_xlim(0, 0.3)
        plt.subplots_adjust(wspace=0, hspace=0)
        ax[0].legend(labels, loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=4, fancybox=True, shadow=True)
        ax[0].text(0.001, 5*1e13, "a", fontsize=30)
        ax[1].text(-0.0015, 1e13, "b", fontsize=30)
        ax[2].text(0.001, 1e9, "c", fontsize=30)
        # plt.show()
        plt.savefig(os.path.join(options.output_folder,'single_island_RT.png'), bbox_inches='tight', pad_inches=0.1)
        plt.close(fig)
    elif action == "plot_jumps":
        full = options.full
        for directory,names in zip(directories, file_names):
            for name in names:
                try:
                    s = SingleResultsProcessor(directory, name, fullOutput=full, vertCurrent=False, graph=False,
                                               reAnalyze=options.re_analyze)
                    if filter(s):
                        print("Plotting file: " + name + " from directory " + directory)
                        fig, fig2, fig3= s.plot_jumps_freq(by_occupation=options.by_occupation,
                                                           by_path_occupation=options.by_path_occupation)
                        fig.canvas.set_window_title(name.replace("_", " "))
                        if options.output_folder:
                            fig.savefig(os.path.join(options.output_folder, name + '_jumps.png'), bbox_inches='tight')
                            plt.close(fig)
                            if fig2 is not None:
                                fig2.savefig(os.path.join(options.output_folder, name + '_jumps_hist.png'), bbox_inches='tight')
                                plt.close(fig2)
                            if fig3 is not None:
                                fig3.savefig(os.path.join(options.output_folder, name + '_jumps_fou.png'),
                                             bbox_inches='tight')
                                plt.close(fig3)
                        else:
                            plt.show()
                except MissingFilesException:
                    print("Missig file for plotting " + name)
                    continue
    elif action == 'plot_score_by_parameter':
        directories_list = []
        for directory, names in zip(directories, file_names):
            directories_list += [directory] * len(names)
        names = [name for files in file_names for name in files]
        m = MultiResultAnalyzer(directories_list, names, out_directory=None, graph=False,reAnalyze=False,
                                relevant_array_params=["Rh", "Rv", "Ch", "Cv"],
                                relevant_running_params=["C_std","R_std", "R_avg","C_avg", "CG_avg", "CG_std", "N", "M", "T"],
                                filter=filter)
        if options.plot3d:
            fig = plt.figure(figsize=FIGSIZE)
            for score in options.scores:
                m.plot_score_by_parameter(score, options.parameters, score, average_results=options.average,
                                          normalizing_parameter_names=options.normalizing, fig=fig)
        else:
            normalizing_idx = 0
            for param in options.parameters:
                fig = plt.figure(figsize=FIGSIZE)
                for score in options.scores:
                    if len(options.normalizing) > 0:
                        m.plot_score_by_parameter(score, [param], score, average_results=options.average,
                                                  normalizing_parameter_names=[options.normalizing[normalizing_idx]],
                                                  fig=fig)
                    else:
                        m.plot_score_by_parameter(score, [param], score, average_results=options.average, fig=fig)
                normalizing_idx += 1
        plt.xlabel(options.xlabel)
        plt.ylabel(options.ylabel)
        plt.legend(["Increasing voltage", "Decreasing voltage"])
        plt.show()
    elif action == 'plot_score_by_disorder':
        directories_list = []
        fig = plt.figure(figsize=FIGSIZE)
        for directory, names in zip(directories, file_names):
            directories_list += [directory] * len(names)
        names = [name for files in file_names for name in files]
        m = MultiResultAnalyzer(directories_list, names, out_directory=None, graph=False, reAnalyze=False,
                                relevant_array_params=["Rh", "Rv", "Ch", "Cv"],
                                relevant_running_params=["C_std","R_std", "R_avg","C_avg", "CG_avg", "CG_std"],
                                filter=filter)
        m.plot_results_by_disorder(options.parameters, options.scores, plot3D=options.plot3d, fig=fig,
                                   average=options.average, window_size=options.window)
        plt.xlabel(options.xlabel)
        plt.ylabel(options.ylabel)
        plt.legend(["Increasing voltage", "Decreasing voltage"])
        plt.show()
    elif action == 'plot_score_by_average':
        directories_list = []
        for directory, names in zip(directories, file_names):
            directories_list += [directory] * len(names)
        names = [name for files in file_names for name in files]
        m = MultiResultAnalyzer(directories_list, names, out_directory=None, graph=False,reAnalyze=options.re_analyze,
                                relevant_array_params=["Rh", "Rv", "Ch", "Cv"],
                                relevant_running_params=["C_std","R_std", "R_avg","C_avg", "CG_avg", "CG_std", "T",
                                                         "M", "N"],
                                filter=filter)
        m.plot_results_by_average(options.parameters, options.scores, plot3D=options.plot3d, average=options.average,
                                  window_size=options.window)
        plt.xlabel(options.xlabel)
        plt.ylabel(options.ylabel)
        #plt.legend(["Increasing voltage", "Decreasing voltage"])
        plt.show()

    elif action == 'plot_scores_together':
        directories_list = []
        for directory, names in zip(directories, file_names):
            directories_list += [directory] * len(names)
        names = [name for files in file_names for name in files]
        m = MultiResultAnalyzer(directories_list, names, out_directory=None, graph=False,reAnalyze=options.re_analyze,
                                relevant_array_params=["Rh", "Rv", "Ch", "Cv"],
                                relevant_running_params=["C_std","R_std", "R_avg","C_avg", "CG_avg", "CG_std"],
                                filter=filter)
        fig, ax = plt.subplots()
        for score in options.scores[1:]:
            m.plot_scores_together(options.scores[0], score, average_results=options.average,
                                   window_size=options.window, fig=fig, ax=ax)
        plt.xlabel(options.xlabel)
        plt.ylabel(options.ylabel)
        if options.plot3d:
            plt.zlabel(options.zlabel)
        plt.show()

    ####### Manual actions ########
    elif action == "compareIV": # compares methods
        output_dir = '/home/kasirershahar/University/Research/Thesis/Figures/'
        directory_graph = "/home/kasirershahar/University/Research/Numerics/jumps_stats"
        directory_gillespie = "/home/kasirershahar/University/Research/simulation_results/single_island/single_island_different_temperature"
        resultNames= ["big_CG_big_R2","big_CG_small_R2","small_CG_big_R2","small_CG_small_R2"]
        namesGraph = ["large_R2_large_CG", "small_R2_large_CG",
                      "large_R2_small_CG", "small_R2_small_CG"]
        namesGillespie = ["array_1_1_big_cg_big_r2_T_", "array_1_1_big_cg_small_r2_T_",
                          "array_1_1_small_cg_big_r2_T_", "array_1_1_small_cg_small_r2_T_"]
        fig, axes = plt.subplots(2,2,figsize=FIGSIZE)
        plot_idx = 0
        for nameGraph, nameGillespie, resultName in zip(namesGraph, namesGillespie, resultNames):
            ax = axes[plot_idx//2, plot_idx%2]
            pGraph = SingleResultsProcessor(directory_graph, nameGraph, fullOutput=False, vertCurrent=False,
                                               reAnalyze=False, graph=True)
            psGillespie = [SingleResultsProcessor(directory_gillespie, nameGillespie + str(T), fullOutput=False, vertCurrent=False,
                                                reAnalyze=True) for T in [0, 0.001, 0.01, 0.1]]
            shift = 0
            psGillespie[0].plot_IV("Dynamical method", err=False, alternative=False, errorevery=1,
                      Vlabel="$V\\frac{C1+C2}{e}$", Ilabel="$I\\frac{\\left(C1+C2\\right)\\left(R1+R2\\right)}{e}$", shift=shift, ax=ax, fig=fig)
            shift += 1
            pGraph.plot_IV("Graph method", err=False, alternative=False,
                           Vlabel="$V\\frac{C1+C2}{e}$", Ilabel="$I\\frac{\\left(C1+C2\\right)\\left(R1+R2\\right)}{e}$", fmt_up='c*', fmt_down='m*', ax=ax, fig=fig)
            for p in psGillespie[1:]:
                p.plot_IV("Dynamical method", err=False, alternative=False, errorevery=1,
                               Vlabel="$V\\frac{C1+C2}{e}$", Ilabel="$I\\frac{\\left(C1+C2\\right)\\left(R1+R2\\right)}{e}$", shift=shift, ax=ax, fig=fig)
                shift += 1
            plot_idx += 1

        axes[0,0].set_xlabel('')
        axes[0, 0].set_xticks([])
        axes[0,1].set_xlabel('')
        axes[0, 1].set_xticks([])
        axes[0, 1].set_ylabel('')
        axes[0, 1].set_yticks([])
        axes[1, 1].set_ylabel('')
        axes[1, 1].set_yticks([])
        plt.subplots_adjust(wspace=0, hspace=0)
        axes[0,0].legend(['Dynamical method (increasing voltage)', 'Dynamical method (decreasing voltage)',
                    'Graph method (increasing voltage)','Graph method (decreasing voltage)'],loc='upper center',
                         bbox_to_anchor=(1, 1.2), ncol=2, fancybox=True, shadow=True)
        axes[0,0].set_zorder(4)
        axes[0,0].text(-0.2,5.8,"a",fontsize=30)
        axes[0,1].set_zorder(3)
        axes[0,1].text(-0.2,5.8,"b",fontsize=30)
        axes[1,0].set_zorder(2)
        axes[1,0].text(-0.2,5.8,"c",fontsize=30)
        axes[1,1].set_zorder(1)
        axes[1,1].text(-0.2,5.8,"d",fontsize=30)
        # plt.show()
        plt.savefig(os.path.join(output_dir, "IV_examples_single_island.png"), bbox_inches = 'tight', pad_inches = 0.01)
        plt.close(fig)

    elif action == "jumps_separation_compare":
        directory = "/home/kasirershahar/University/Research/Numerics/jumps_stats"
        plt.figure(figsize=FIGSIZE)
        shift = 0
        for Cratio in [0.1,0.2,0.3,0.4,0.5,0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
            Rratios = [0.1,0.2,0.3,0.4,0.5,3,4,5,6,7,8,9,10]
            names = ["single_island_Cratio_" + str(Cratio) + "_Rratio_" + str(Rratio) for Rratio in Rratios]
            full = True
            save_re_analysis = False
            directories_list = [directory]*len(names)
            m = MultiResultAnalyzer(directories_list, names, out_directory=None, graph=True)
            labelprefix = "C2/C1 = " + str(Cratio)
            Rratios = np.array(Rratios)
            CG=2
            C1 = 1/(Cratio + 1)
            C2 = Cratio/(Cratio + 1)
            approx = []
            for Rratio in Rratios:
                R1 = 1/(Rratio + 1)
                R2 = Rratio/(Rratio + 1)
                approx.append(approx_jumps_separation(R1, R2, C1, C2, CG))
            m.plot_score("jumpSeparationUp", xs=Rratios, xlabel="R2/R1", ylabel="Average voltage difference [e/(C1+C2)]", fmt='r.',
                          shift=shift)
            m.plot_score("jumpSeparationDown", xs=Rratios, xlabel="R2/R1", ylabel="Average voltage difference [e/(C1+C2)]", fmt='b.',
                          shift=shift)
            plt.plot(Rratios, np.array(approx) + shift, 'orange')
            shift += 1
        plt.show()

    elif action == "jumps_up_down_compare":
        directory = "/home/kasirershahar/University/Research/Numerics/jumps_stats"
        plt.figure(figsize=FIGSIZE)
        shift = 0
        for Rratio in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,2,3,4,5,6,7,8,9,10]:
            Cratios = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10]
            names = ["single_island_Cratio_" + str(Cratio) + "_Rratio_" + str(Rratio) for Cratio in Cratios]
            full = True
            save_re_analysis = False
            directories_list = [directory]*len(names)
            m = MultiResultAnalyzer(directories_list, names, out_directory=None, graph=True)
            loopV, _, _ = m.get_relevant_scores('firstJumpV')
            score,_,_ = m.get_relevant_scores("jumpSeparationUpDown")
            loopV = np.array(loopV)
            loopV[score==0] = 0
            labelprefix = "R2 < R1" if Rratio == 0.1 else "R2 > R1"
            Cratios = np.array(Cratios)
            CG=2
            VG=0.5
            approx = []
            for i,Cratio in enumerate(Cratios):
                C1 = 1/(Cratio + 1)
                C2 = Cratio/(Cratio + 1)
                R1 = 1/(Rratio + 1)
                R2 = Rratio/(Rratio + 1)
                approx.append(approx_jump_up_down(R1, R2, C1, C2, CG, loopV[i]))

            m.plot_score("jumpSeparationUpDown", xs=Cratios, xlabel="C2/C1", ylabel="Average voltage difference [e/(C1+C2)]",
                         label=labelprefix + " numerical results", shift=shift, fmt='m.')
            plt.plot(Cratios[loopV > 0], np.array(approx)[loopV > 0]+shift, "orange", label=labelprefix + " analytic approximation")
            shift += 1
        plt.show()

    elif action == "jumps_height_compare":
        directory = "/home/kasirershahar/University/Research/Numerics/jumps_stats"
        plt.figure(figsize=FIGSIZE)
        shift = 0
        for Rratio in [0.1, 0.2, 0.3, 0.4, 2,3,4,5,6,7,8,9,10]:
            Cratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            names = ["single_island_Cratio_" + str(Cratio) + "_Rratio_" + str(Rratio) for Cratio in Cratios]
            full = True
            save_re_analysis = False
            directories_list = [directory] * len(names)
            m = MultiResultAnalyzer(directories_list, names, out_directory=None, graph=True)
            labelprefix = "R2 < R1" if Rratio == 0.1 else "R2 > R1"
            Cratios = np.array(Cratios)
            CG = 2
            C1 = 1 / (Cratios + 1)
            C2 = Cratios / (Cratios + 1)
            R1 = 1 / (Rratio + 1)
            R2 = Rratio / (Rratio + 1)
            approx = 1/((C1+C2+CG)*(R1+R2)) + shift
            m.plot_score("jumpHeightUp", xs=Cratios, xlabel="C2/C1",
                         ylabel="Average current difference [e/(C1+C2)(R1+R2)]",
                         label=labelprefix + " increasing voltage numerical results", shift=shift, fmt=".r")
            m.plot_score("jumpHeightDown", xs=Cratios, xlabel="C2/C1",
                         ylabel="Average current difference [e/(C1+C2)(R1+R2)]",
                         label=labelprefix + " decreasing voltage numerical results", shift=shift, fmt=".b")
            plt.plot(Cratios, approx, "orange",label=labelprefix + " analytic approximation")
            shift += 1
        plt.show()

    elif action == "hysteresis_area_compare":
        directory = "/home/kasirershahar/University/Research/Numerics/jumps_stats"
        plt.figure(figsize=FIGSIZE)
        shift = 0
        for Cratio in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2, 3, 4, 5, 6, 7, 8, 9, 10]:
            Rratios = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2, 3, 4, 5, 6, 7, 8, 9, 10]
            names = ["single_island_Cratio_" + str(Cratio) + "_Rratio_" + str(Rratio) for Rratio in Rratios]
            full = True
            save_re_analysis = False
            directories_list = [directory] * len(names)
            m = MultiResultAnalyzer(directories_list, names, out_directory=None, graph=True)
            CG = 2
            VG = 0.5
            approx = []

            for Rratio in Rratios:
                C1 = 1 / (Cratio + 1)
                C2 = Cratio / (Cratio + 1)
                R1 = 1 / (Rratio + 1)
                R2 = Rratio / (Rratio + 1)
                approx.append(approx_hysteresis_total_area(R1, R2, C1, C2, CG, VG, 4))
            m.plot_score("hysteresisArea", xs=Rratios, xlabel="R2/R1",
                         ylabel="Average loop area [e^2/(C1+C2)^2(R1+R2)]", shift=shift, fmt='m.')
            plt.plot(Rratios, np.array(approx) + shift, 'orange', linestyle='dashed')
            shift += 1
        plt.show()

    elif action == "threshold_voltage_compare":
        directory = "/home/kasirershahar/University/Research/Numerics/jumps_stats"
        plt.figure(figsize=FIGSIZE)
        shift = 0
        for Cratio in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
            Rratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            names = ["single_island_Cratio_" + str(Cratio) + "_Rratio_" + str(Rratio) for Rratio in Rratios]
            full = True
            save_re_analysis = False
            directories_list = [directory] * len(names)
            m = MultiResultAnalyzer(directories_list, names, out_directory=None, graph=True)
            labelprefix = "C2 = C1" if Cratio == 1 else "C2 = 2C1"
            approx = []
            CG = 2
            VG = 0.5
            nini = 0
            for Rratio in Rratios:
                C1 = 1 / (Cratio + 1)
                C2 = Cratio / (Cratio + 1)
                R1 = 1 / (Rratio + 1)
                R2 = Rratio / (Rratio + 1)
                approx.append(approx_threshold_voltage(R1, R2, C1, C2, CG, VG, nini))
            m.plot_score("thresholdVoltageUp", xs=Rratios, xlabel="R2/R1",
                         ylabel="Threshold voltage [e^2/(C1+C2)^2(R1+R2)]",
                         label=labelprefix + " numerical results", shift=shift, fmt=".m")
            plt.plot(Rratios, np.array(approx)+shift, "orange", label=labelprefix + " analytic approximation")
            shift+=1
        plt.show()

    elif action == "jumps_num_compare":
        directory = "/home/kasirershahar/University/Research/Numerics/jumps_stats"
        plt.figure(figsize=FIGSIZE)
        Rratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        Cratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6]
        names = ["single_island_Cratio_" + str(Cratio) + "_Rratio_" + str(Rratio) for Cratio in Cratios for Rratio in Rratios]
        full = True
        save_re_analysis = False
        directories_list = [directory] * len(names)
        m = MultiResultAnalyzer(directories_list, names, out_directory=None, graph=True)
        shift = 0
        CG = 2
        VG = 0.5
        nini = 0
        groups = np.arange(len(names))
        m.plot_score_by("jumpNumUp", xs=np.array(Cratios), xlabel="C2/C1",
                     ylabel="Number of jumps",
                     label="numerical results", shift=shift)
        m.plot_score("jumpNumDown", xs=np.array(Cratios), xlabel="C2/C1",
                     ylabel="Number of jumps",
                     label="numerical results", shift=shift)
        # plt.plot(Rratios, np.array(approx) + shift, "orange", label=labelprefix + " analytic approximation")
        plt.show()

    elif action == 'IV_different_temperatures':
        directory = "/home/kasirershahar/University/Research/old_results/same_array_different_temperature"
        Ts = [0.001, 0.002, 0.005, 0.007, 0.01, 0.015, 0.02]
        names = ["array_10_10_T_" + str(T) for T in Ts]
        shift=0
        for T,name in zip(Ts,names):
            p = SingleResultsProcessor(directory, name, reAnalyze=True, graph=False, fullOutput=True)
            p.plot_IV("T = " + str(T), Vnorm=1/2, Inorm=1/20, shift=shift, err=True, errorevery=5,alternative=False,
                      Ilabel="I(<R><C>)/e",Vlabel="V<C>/e")
            shift+=0.3
        plt.xlim(1,2)
        plt.ylim(-0.1,2.5)
        plt.legend(["Increasing voltage", "Decreasing voltage"])
        plt.show()

    elif action == 'score_by_cg':
        directories_list = []
        for directory, names in zip(directories, file_names):
            directories_list += [directory] * len(names)
        names = [name for files in file_names for name in files]
        m = MultiResultAnalyzer(directories_list, names, out_directory=None, graph=False,reAnalyze=options.re_analyze,
                                relevant_array_params=["Rh", "Rv", "Ch", "Cv"],
                                relevant_running_params=["C_std","R_std", "R_avg","C_avg", "CG_avg", "CG_std", "T",
                                                         "M", "N"],
                                filter=filter)

        fig, ax1 = plt.subplots()
        m.plot_results_by_average(["CG"], ["jumpSeparationUp"], plot3D=options.plot3d, average=options.average,
                                  window_size=options.window, fmt='.b', fig=fig, ax=ax1)
        m.plot_results_by_average(["CG"], ["jumpSeparationDown"], plot3D=options.plot3d, average=options.average,
                                  window_size=options.window, fmt='*b', fig=fig, ax=ax1)
        ax1.set_xlabel("$\\frac{C_G}{\\left<C\\right>}$")
        ax1.set_ylabel("$\\Delta V \\frac{\\left<C\\right>}{e}$", color='blue')
        ax1.tick_params(axis='y', labelcolor='blue')
        ax2 = ax1.twinx()
        m.plot_results_by_average(["CG"], ["jumpHeightUp"], plot3D=options.plot3d, average=options.average,
                                  window_size=options.window, fmt='.r', fig=fig, ax=ax2)
        m.plot_results_by_average(["CG"], ["jumpHeightDown"], plot3D=options.plot3d, average=options.average,
                                  window_size=options.window, fmt='*r', fig=fig, ax=ax2)
        ax2.set_ylabel("$\\Delta I \\frac{\\left<C\\right>\\left<R\\right>}{e}$", color='red')
        ax2.tick_params(axis='y', labelcolor='red')
        ax3 = ax1.twinx()
        m.plot_results_by_average(["CG"], ["hysteresisArea"], plot3D=options.plot3d, average=options.average,
                                  window_size=options.window, fmt='om', fig=fig, ax=ax3)
        ax3.set_ylabel("Hysteresis area $\\left[\\frac{e^2}{\\left<C\\right>^2\\left<R\\right>}\\right]$", color='m')
        ax3.tick_params(axis='y', labelcolor='m')
        ax3.spines['right'].set_position(('outward', 90))
        fig.tight_layout()
        plt.show()

    elif action == 'score_by_disorder_zero_T':
        c_window=0.05
        r_window=0.2
        directories_list = []
        for directory, names in zip(directories, file_names):
            directories_list += [directory] * len(names)
        names = [name for files in file_names for name in files]
        m = MultiResultAnalyzer([directories[1]]*len(file_names[1]), file_names[1], out_directory=None, graph=False,reAnalyze=options.re_analyze,
                                relevant_array_params=["Rh", "Rv", "Ch", "Cv"],
                                relevant_running_params=["C_std","R_std", "R_avg","C_avg", "CG_avg", "CG_std", "T",
                                                         "M", "N"],
                                filter=filter)

        fig, axes = plt.subplots(3,2,figsize=FIGSIZE)
        ax1 = axes[0,0]
        m.plot_results_by_disorder(["C"], ["jumpSeparationUp"], plot3D=options.plot3d, average=True,
                                  window_size=c_window, fmt='^b', fig=fig, ax=ax1)
        m.plot_results_by_disorder(["C"], ["jumpSeparationDown"], plot3D=options.plot3d, average=True,
                                  window_size=c_window, fmt='vb', fig=fig, ax=ax1)

        ax1.set_ylabel("$\\Delta V \\frac{\\left<C\\right>}{e}$", color='blue')
        ax1.tick_params(axis='y', labelcolor='blue')
        ax2 = axes[1,0]
        m.plot_results_by_disorder(["C"], ["jumpHeightUp"], plot3D=options.plot3d, average=True,
                                  window_size=c_window, fmt='^r', fig=fig, ax=ax2)
        m.plot_results_by_disorder(["C"], ["jumpHeightDown"], plot3D=options.plot3d, average=True,
                                  window_size=c_window, fmt='vr', fig=fig, ax=ax2)
        ax2.set_ylabel("$\\Delta I \\frac{\\left<C\\right>\\left<R\\right>}{e}$", color='red')
        ax2.tick_params(axis='y', labelcolor='red')
        ax3 = axes[2,0]
        m.plot_results_by_disorder(["C"], ["hysteresisArea"], plot3D=options.plot3d, average=True,
                                  window_size=c_window, fmt='om', fig=fig, ax=ax3)
        ax3.set_ylabel("Hysteresis area $\\left[\\frac{e^2}{\\left<C\\right>^2\\left<R\\right>}\\right]$", color='m')
        ax3.tick_params(axis='y', labelcolor='m')
        ax3.set_xlabel("$\\sigma_C/\\left<C\\right>$")
        ax1.set_xticklabels([])
        ax2.set_xticklabels([])
        m = MultiResultAnalyzer(directories_list, names, out_directory=None, graph=False,
                                reAnalyze=options.re_analyze,
                                relevant_array_params=["Rh", "Rv", "Ch", "Cv"],
                                relevant_running_params=["C_std", "R_std", "R_avg", "C_avg", "CG_avg", "CG_std", "T",
                                                         "M", "N"],
                                filter=filter)
        ax4 = axes[0,1]
        m.plot_results_by_disorder(["R"], ["jumpSeparationUp"], plot3D=options.plot3d, average=True,
                                   window_size=r_window, fmt='^b', fig=fig, ax=ax4)
        m.plot_results_by_disorder(["R"], ["jumpSeparationDown"], plot3D=options.plot3d, average=True,
                                   window_size=r_window, fmt='vb', fig=fig, ax=ax4)


        ax4.tick_params(axis='y', labelcolor='blue')
        ax5 = axes[1,1]
        m.plot_results_by_disorder(["R"], ["jumpHeightUp"], plot3D=options.plot3d, average=True,
                                   window_size=r_window, fmt='^r', fig=fig, ax=ax5)
        m.plot_results_by_disorder(["R"], ["jumpHeightDown"], plot3D=options.plot3d, average=True,
                                   window_size=r_window, fmt='vr', fig=fig, ax=ax5)
        ax5.tick_params(axis='y', labelcolor='red')
        ax6 = axes[2,1]
        m.plot_results_by_disorder(["R"], ["hysteresisArea"], plot3D=options.plot3d, average=True,
                                   window_size=r_window, fmt='om', fig=fig, ax=ax6)
        ax6.tick_params(axis='y', labelcolor='m')
        ax6.set_xlabel("$\\sigma_R/\\left<R\\right>$")
        ax4.set_xticklabels([])
        ax5.set_xticklabels([])
        axes[0, 0].text(0.04, 0.21, "a", fontsize=30)
        axes[1, 0].text(0.04, 0.0061, "b", fontsize=30)
        axes[2, 0].text(0.04, 0.0012, "c", fontsize=30)
        axes[0, 1].text(0.08, 0.215, "d", fontsize=30)
        axes[1, 1].text(0.08, 0.53, "e", fontsize=30)
        axes[2, 1].text(0.08, 0.255, "f", fontsize=30)
        if options.output_folder:
            fig.savefig(os.path.join(options.output_folder,  'hysteresis_and_jumps_by_disorders_zero_T.png'), bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()

    elif action == 'score_by_r_disorder':
        zero_temp_window=0.2
        finite_temp_window=0.1
        m = MultiResultAnalyzer([directories[0]]*len(file_names[0]) + [directories[1]]*len(file_names[1]) ,
                                file_names[0] + file_names[1],
                                out_directory=None, graph=False,reAnalyze=options.re_analyze,
                                relevant_array_params=["Rh", "Rv", "Ch", "Cv"],
                                relevant_running_params=["C_std","R_std", "R_avg","C_avg", "CG_avg", "CG_std", "T",
                                                         "M", "N"],
                                filter=filter)

        fig, axes = plt.subplots(4,2,figsize=FIGSIZE)
        ax1 = axes[0,0]
        m.plot_results_by_disorder(["R"], ["jumpSeparationUp"], plot3D=options.plot3d, average=True,
                                  window_size=zero_temp_window, fmt='^b', fig=fig, ax=ax1)
        m.plot_results_by_disorder(["R"], ["jumpSeparationDown"], plot3D=options.plot3d, average=True,
                                  window_size=zero_temp_window, fmt='vb', fig=fig, ax=ax1)

        ax1.set_ylabel("$\\Delta V \\frac{\\left<C\\right>}{e}$", color='blue')
        ax1.tick_params(axis='y', labelcolor='blue')
        ax2 = axes[1,0]
        m.plot_results_by_disorder(["R"], ["jumpHeightUp"], plot3D=options.plot3d, average=True,
                                  window_size=zero_temp_window, fmt='^r', fig=fig, ax=ax2)
        m.plot_results_by_disorder(["R"], ["jumpHeightDown"], plot3D=options.plot3d, average=True,
                                  window_size=zero_temp_window, fmt='vr', fig=fig, ax=ax2)
        ax2.set_ylabel("$\\Delta I \\frac{\\left<C\\right>\\left<R\\right>}{e}$", color='red')
        ax2.tick_params(axis='y', labelcolor='red')
        ax3 = axes[2,0]
        m.plot_results_by_disorder(["R"], ["hysteresisArea"], plot3D=options.plot3d, average=True,
                                  window_size=zero_temp_window, fmt='om', fig=fig, ax=ax3)
        ax3.set_ylabel("Hysteresis area $\\left[\\frac{e^2}{\\left<C\\right>^2\\left<R\\right>}\\right]$", color='m')
        ax3.tick_params(axis='y', labelcolor='m')
        ax3.set_xlabel("$\\sigma_R/\\left<R\\right>$")
        ax4 = axes[3, 0]
        m.plot_results_by_disorder(["R"], ["thresholdVoltageUp"], plot3D=options.plot3d, average=True,
                                   window_size=zero_temp_window, fmt='^g', fig=fig, ax=ax4)
        m.plot_results_by_disorder(["R"], ["thresholdVoltageDown"], plot3D=options.plot3d, average=True,
                                   window_size=zero_temp_window, fmt='vg', fig=fig, ax=ax4)
        ax4.set_ylabel("$V^{th} \\frac{\\left<C\\right>}{e}$", color='g')
        ax4.tick_params(axis='y', labelcolor='g')
        ax4.set_xlabel("$\\sigma_R/\\left<R\\right>$")
        ax1.set_xticklabels([])
        ax2.set_xticklabels([])
        ax3.set_xticklabels([])
        m = MultiResultAnalyzer([directories[2]]*len(file_names[2]), file_names[2], out_directory=None, graph=False,
                                reAnalyze=options.re_analyze,
                                relevant_array_params=["Rh", "Rv", "Ch", "Cv"],
                                relevant_running_params=["C_std", "R_std", "R_avg", "C_avg", "CG_avg", "CG_std", "T",
                                                         "M", "N"],
                                filter=filter)
        ax5 = axes[0,1]
        m.plot_results_by_disorder(["R"], ["jumpSeparationUp"], plot3D=options.plot3d, average=True,
                                   window_size=finite_temp_window, fmt='^b', fig=fig, ax=ax5)
        m.plot_results_by_disorder(["R"], ["jumpSeparationDown"], plot3D=options.plot3d, average=True,
                                   window_size=finite_temp_window, fmt='vb', fig=fig, ax=ax5)


        ax5.tick_params(axis='y', labelcolor='blue')
        ax6 = axes[1,1]
        m.plot_results_by_disorder(["R"], ["jumpHeightUp"], plot3D=options.plot3d, average=True,
                                   window_size=finite_temp_window, fmt='^r', fig=fig, ax=ax6)
        m.plot_results_by_disorder(["R"], ["jumpHeightDown"], plot3D=options.plot3d, average=True,
                                   window_size=finite_temp_window, fmt='vr', fig=fig, ax=ax6)
        ax6.tick_params(axis='y', labelcolor='red')
        ax7 = axes[2,1]
        m.plot_results_by_disorder(["R"], ["hysteresisArea"], plot3D=options.plot3d, average=True,
                                   window_size=finite_temp_window, fmt='om', fig=fig, ax=ax7)
        ax7.tick_params(axis='y', labelcolor='m')
        ax7.set_xlabel("$\\sigma_R/\\left<R\\right>$")
        ax8 = axes[3, 1]
        m.plot_results_by_disorder(["R"], ["thresholdVoltageUp"], plot3D=options.plot3d, average=True,
                                   window_size=finite_temp_window, fmt='^g', fig=fig, ax=ax8)
        m.plot_results_by_disorder(["R"], ["thresholdVoltageDown"], plot3D=options.plot3d, average=True,
                                   window_size=finite_temp_window, fmt='vg', fig=fig, ax=ax8)
        ax8.set_ylabel("$V^{th} \\frac{\\left<C\\right>}{e}$", color='g')
        ax8.tick_params(axis='y', labelcolor='g')
        ax8.set_xlabel("$\\sigma_R/\\left<R\\right>$")
        ax5.set_xticklabels([])
        ax6.set_xticklabels([])
        ax7.set_xticklabels([])
        xmin, xmax, _, _ = axes[1, 1].axis()
        axes[0, 1].set_xlim(xmin, xmax)
        plt.subplots_adjust(wspace=0.13, hspace=0)
        add_text_upper_left_corner(axes[0, 0], "a")
        add_text_upper_left_corner(axes[1, 0], "b")
        add_text_upper_left_corner(axes[2, 0], "c")
        add_text_upper_left_corner(axes[3, 0], "d")
        add_text_upper_left_corner(axes[0, 1], "e")
        add_text_upper_left_corner(axes[1, 1], "f")
        add_text_upper_left_corner(axes[2, 1], "g")
        add_text_upper_left_corner(axes[3, 1], "h")
        if options.output_folder:
            fig.savefig(os.path.join(options.output_folder,  'hysteresis_and_jumps_by_r_disorders.png'), bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()

    elif action == 'IV_path_compare':
        vertical_directory = '/home/kasirershahar/University/Research/simulation_results/zero_temperature/bgu_1d_vertical_arrays/'
        single_directory = '/home/kasirershahar/University/Research/simulation_results/single_island/single_island_for_vertical_arrays/'
        vertical_files = ['array_2_1_custome_r_c_disorder_cg_1_run_1', 'array_2_1_custome_r_c_disorder_cg_5_run_1',
                          'array_3_1_custome_r_c_disorder_cg_10_run_1', 'array_3_1_custome_r_c_disorder_cg_50_run_1']
        single_files = [['array_2_1_cg_1_row_0','array_2_1_cg_1_row_1'],['array_2_1_cg_5_row_0','array_2_1_cg_5_row_1'],
                        ['array_3_1_cg_10_row_0','array_3_1_cg_10_row_1','array_3_1_cg_10_row_2'],
                        ['array_3_1_cg_50_row_0','array_3_1_cg_50_row_1','array_3_1_cg_50_row_2']]
        fmt_list = ['r--','m--','r--','m--','r--','g--','m--','r--','g--','m--']
        fmt_idx = 0
        cg_list = [1,5,10,50]
        cg_idx = 0
        ax_idx = 0
        fig, axes = plt.subplots(2, 2, figsize=FIGSIZE)
        for vertical_file,single_files_list in zip(vertical_files, single_files):
            ax = axes[ax_idx//2,ax_idx%2]
            p = SingleResultsProcessor(vertical_directory, vertical_file, reAnalyze=False, graph=False, fullOutput=True)
            p.plot_current_by_paths(fig, ax)
            p.createCapacitanceMatrix()
            vert_invC = p.invC
            for row, file in enumerate(single_files_list):
                ps = SingleResultsProcessor(single_directory, file, reAnalyze=False, graph=False, fullOutput=True)
                ps.createCapacitanceMatrix()
                single_invC = ps.invC
                c_factor = 1/(2*vert_invC[row, row])
                ps.plot_IV(ax, fig, Inorm=c_factor, Vnorm=c_factor, fmt_up=fmt_list[fmt_idx], fmt_down=fmt_list[fmt_idx])
                fmt_idx += 1

            if ax_idx < 2:
                ax.set_xlabel("")
                ax.set_xticklabels([])
            else:
                ax.set_xticklabels([0,0.25,0.5,0.75,1,1.25,1.5,1.75])
                ax.set_xlabel("$V\\frac{\\left<C\\right>}{e}$")
            if ax_idx%2 > 0:
                ax.set_ylabel("")
                ax.set_yticklabels([])
            else:
                ax.set_ylabel("$I\\frac{\\left<R\\right>\\left<C\\right>}{e}$")
            ax.set_title("$\\frac{C_G}{\\left<C\\right>} = %d$" % cg_list[cg_idx],position=(0.5, 0.8))
            ax.set_xlim(0,2)
            ax.set_ylim(-0.05,1)
            cg_idx += 1
            ax_idx += 1
        plt.subplots_adjust(wspace=0, hspace=0)
        axes[0, 0].text(0.01, 0.92, "a", fontsize=50)
        axes[0, 1].text(0.01, 0.92, "b", fontsize=50)
        axes[1, 0].text(0.01, 0.92, "c", fontsize=50)
        axes[1, 1].text(0.01, 0.92, "d", fontsize=50)
        fig.savefig(os.path.join(options.output_folder,  'vertical_single_compare.png'), bbox_inches='tight')

    elif action == 'plot_scores_with_threshold_voltage':
        fig, ax = plt.subplots(2, 2, figsize=FIGSIZE)
        directory_index=0
        windows = [0.02, 0.05]
        for directory in directories:
            names = file_names[directory_index]
            directories_list = [directory] *len(names)
            m = MultiResultAnalyzer(directories_list, names, out_directory=None, graph=False,reAnalyze=options.re_analyze,
                                    relevant_array_params=["Rh", "Rv", "Ch", "Cv"],
                                    relevant_running_params=["C_std","R_std", "R_avg","C_avg", "CG_avg", "CG_std"],
                                    filter=filter)

            m.plot_scores_together(options.scores[0], options.scores[1], average_results=options.average,
                                       window_size=windows[directory_index], fig=fig, ax=ax[0,directory_index],fmt='.b')
            m.plot_scores_together(options.scores[0], options.scores[2], average_results=options.average,
                                     window_size=windows[directory_index], fig=fig, ax=ax[1,directory_index],fmt='.m')
            directory_index += 1
        ax[1,0].set_xlabel("$V^{th}\\frac{\\left<C\\right>}{e}$")
        ax[1, 1].set_xlabel("$V^{th}\\frac{\\left<C\\right>}{e}$")
        ax[0,0].set_ylabel("$\\Delta I \\frac{\\left<C\\right>\\left<R\\right>}{e}$", color='b')
        ax[1, 0].set_ylabel("Hysteresis area $\\left[\\frac{e^2}{\\left<C\\right>^2\\left<R\\right>}\\right]$", color='m')
        ax[0,0].tick_params(axis='y', labelcolor='b')
        ax[0, 1].tick_params(axis='y', labelcolor='b')
        ax[1,0].tick_params(axis='y', labelcolor='m')
        ax[1, 1].tick_params(axis='y', labelcolor='m')
        ax[0,0].set_xticklabels([])
        ax[0, 1].set_xticklabels([])
        ax[0, 0].text(0.88, 0.025, "a", fontsize=30)
        ax[1, 0].text(0.88, 0.00052, "b", fontsize=30)
        ax[0, 1].text(0.415, 0.66, "c", fontsize=30)
        ax[1, 1].text(0.415, 0.42, "d", fontsize=30)
        plt.subplots_adjust(wspace=0.13, hspace=0)
        if options.output_folder:
            fig.savefig(os.path.join(options.output_folder, 'jumps_and_hysteresis_by_threshold.png'), bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()

    elif action == 'plot_small_v_fit':
        for directory, names in zip(directories, file_names):
            full = True
            for name in names:
                s = SingleResultsProcessor(directory, name, fullOutput=full, reAnalyze=False)
                s.power_law_fit(err=True, errorevery=10, plot=True)
                plt.show()
    elif action == 'plot_current_by_path':
        for directory, names in zip(directories, file_names):
            for name in names:
                fig = plt.figure(figsize=FIGSIZE)
                s = SingleResultsProcessor(directory, name, vertCurrent=False, reAnalyze=False)
                s.plot_current_by_paths(fig=fig)
                plt.show()
    elif action == 'plot_jumps_for_2d_vertical_arrays':
        fig, axes = plt.subplots(2,2,figsize=FIGSIZE)
        plt.subplots_adjust(wspace=0, hspace=0)
        directory = '/home/kasirershahar/University/Research/simulation_results/zero_temperature/bgu_1d_vertical_arrays/'
        names = ['array_2_1_custome_r_c_disorder_cg_1_run_1','array_2_1_custome_r_c_disorder_cg_5_run_1',
                 'array_3_1_custome_r_c_disorder_cg_10_run_1','array_3_1_custome_r_c_disorder_cg_50_run_1']
        for index,name in enumerate(names):
            ax = axes[index//2,index%2]
            s = SingleResultsProcessor(directory, name, fullOutput=True, vertCurrent=False, graph=False,
                                       reAnalyze=options.re_analyze)
            fig, fig2, fig3 = s.plot_jumps_freq(by_occupation=options.by_occupation,
                              by_path_occupation=options.by_path_occupation, ax1=ax, fig=fig)
            plt.close(fig2)
            plt.close(fig3)
        axes[0, 0].set_xlabel("")
        axes[0, 0].set_xticklabels([])
        axes[0, 1].set_xlabel("")
        axes[0, 1].set_xticklabels([])
        axes[0, 1].set_ylabel("")
        axes[0, 1].set_yticklabels([])
        axes[1, 1].set_ylabel("")
        axes[1, 1].set_yticklabels([])
        axes[0, 0].text(-0.08, 1.5, "a", fontsize=30)
        axes[0, 1].text(-0.08, 1.5, "b", fontsize=30)
        axes[1, 0].text(-0.08, 2.3, "c", fontsize=30)
        axes[1, 1].text(-0.08, 2.3, "d", fontsize=30)
        if options.output_folder:
            fig.savefig(os.path.join(options.output_folder, 'jumps_1d_vertical_array.png'), bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()
    elif action == 'plot_jumps_for_resistance_effect':
        fig, axes = plt.subplots(3,2,figsize=FIGSIZE)
        plt.subplots_adjust(wspace=0, hspace=0)
        directories = ['/home/kasirershahar/University/Research/simulation_results/zero_temperature/bgu_1d_horizontal_arrays']*3+\
                      ['/home/kasirershahar/University/Research/simulation_results/zero_temperature/bgu_jumps_check']*3
        names = ['array_1_2_no_disorder_cg_5', 'array_1_2_lambda_shape_cg_10', 'array_1_2_u_shape_cg_10',
                 'array_5_5_const_r_c_disorder_cg_10_run_3', 'array_5_5_const_r_c_disorder_cg_10_run_3_v1',
                 'array_5_5_const_r_c_disorder_cg_10_run_3_v3']
        index=0
        for directory,name in zip(directories,names):
            ax = axes[index%3,index//3]
            s = SingleResultsProcessor(directory, name, fullOutput=True, vertCurrent=False, graph=False,
                                       reAnalyze=options.re_analyze)

            fig, fig2, fig3 = s.plot_jumps_freq(by_occupation=index<3,
                              by_path_occupation=index >= 3, ax1=ax, fig=fig)
            plt.close(fig2)
            plt.close(fig3)
            index += 1
        axes[0, 0].set_xlabel("")
        axes[0, 0].set_xticklabels([])
        axes[0, 1].set_xlabel("")
        axes[0, 1].set_xticklabels([])
        axes[1, 0].set_xlabel("")
        axes[1, 0].set_xticklabels([])
        axes[1, 1].set_xlabel("")
        axes[1, 1].set_xticklabels([])
        axes[0, 1].set_ylabel("")
        axes[0, 1].set_yticklabels([])
        axes[1, 1].set_ylabel("")
        axes[1, 1].set_yticklabels([])
        axes[2, 1].set_ylabel("")
        axes[2, 1].set_yticklabels([])
        axes[0, 0].text(-0.08, 0.35, "a", fontsize=30)
        axes[1, 0].text(-0.08, 0.35, "b", fontsize=30)
        axes[2, 0].text(-0.08, 0.35, "c", fontsize=30)
        axes[0, 1].text(0.32, 0.35, "d", fontsize=30)
        axes[1, 1].text(0.32, 0.35, "e", fontsize=30)
        axes[2, 1].text(0.32, 0.35, "f", fontsize=30)
        for ax in axes.flatten():
            ax.set_ylim(-0.01,0.39)
        for ax in axes[:,1]:
            ax.set_xlim(0.3,3.51)
        if options.output_folder:
            fig.savefig(os.path.join(options.output_folder, 'resistance_effect.png'), bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()
    elif action == 'plot_results':
         full = True
         for directory, names in zip(directories, file_names):
             for name in names:
                 s = SingleResultsProcessor(directory, name, fullOutput=full, vertCurrent=False, reAnalyze=True)
                 # s.plot_jumps_freq("I")
                 s.plot_results()
                 s.plot_array_params("R")
                 s.plot_array_params("C")
                 # s.plot_voltage()
                 # s.plot_IV("IV",err=False, errorevery=10, alternative=True)
                 plt.show()

    elif action == 'print_scores':
         full = True
         for directory, names in zip(directories, file_names):
             for name in names:
                 s = SingleResultsProcessor(directory, name, fullOutput=full, vertCurrent=False, reAnalyze=True)
                 if filter(s):
                    s.print_scores(for_paths=True)

    else:
        print("Unknown action")
        exit(0)
