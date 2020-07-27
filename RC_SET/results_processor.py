import numpy as np
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
font = {'family' : 'sans-serif',
        'weight' : 'bold',
        'size' : 22}
matplotlib.rc('font', **font)
from matplotlib import colors
from scipy import integrate
from ast import literal_eval
import os, sys
import re
from sklearn.mixture import GaussianMixture
from mpl_toolkits.mplot3d import Axes3D

EPS = 1e-6
FIGSIZE=(30,16)

def flattenToColumn(a):
    return a.reshape((a.size, 1))

def flattenToRow(a):
    return a.reshape((1, a.size))

class SingleResultsProcessor:
    """ Used for proccessing result from one dot array simulation"""
    def __init__(self, directory, file_name, fullOutput=False, vertCurrent=False):
        self.basePath = os.path.join(directory, file_name)
        self.fileName = file_name
        self.directory = directory
        self.full = fullOutput
        self.vert = vertCurrent
        self.I = None
        self.V = None
        self.IErr = None
        self.n = None
        self.Q = None
        self.nErr = None
        self.QErr = None
        self.mid_idx = 0
        self.runningParams = dict()
        self.arrayParams = dict()
        self.load_results()
        self.load_params()
        self.rows = self.runningParams["M"]
        self.columns = self.runningParams["N"]

    def reAnalyzeI(self):
        full_I = self.full_I.T
        self.alternativeI = []
        self.alternativeIErr = []
        self.alternativeV = []
        for i in range(len(full_I)):
            min_idx = max(0, i-5)
            max_idx = min(len(full_I)-1, i+5)
            point = np.hstack([full_I[j] for j in range(min_idx, max_idx)])
            bins = bistabilityAnalysis(point)
            # if i > 450:
            #     plt.hist(point)
            #     plt.show()
            if i > 0 and len(bins) > 1 and len(bins[0]) > 0 and len(bins[1]) > 0:
                avg1 = np.average(bins[0])
                avg2 = np.average(bins[1])
                if np.abs(avg1 - self.I[i-1]) < np.abs(avg2 - self.I[i-1]):
                # if len(bins[0]) > len(bins[1]):
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
        return True


    def load_results(self):
        I_file = os.path.join(self.basePath + '_I.npy')
        V_file = os.path.join(self.basePath + '_V.npy')
        IErr_file = os.path.join(self.basePath + '_IErr.npy')
        self.I = np.load(I_file)
        self.IErr = np.load(IErr_file)
        self.V = np.load(V_file)
        self.mid_idx = self.V.size // 2
        if self.full:
            n_file = os.path.join(self.basePath + '_n.npy')
            Q_file = os.path.join(self.basePath + '_Q.npy')
            nErr_file = os.path.join(self.basePath + '_nErr.npy')
            QErr_file = os.path.join(self.basePath + '_QErr.npy')
            full_I_file = os.path.join(self.basePath + '_full_I.npy')
            self.n = np.load(n_file)
            self.Q = np.load(Q_file)
            self.nErr = np.load(nErr_file)
            self.QErr = np.load(QErr_file)
            self.full_I = np.load(full_I_file)
            self.reAnalyzeI()
        if self.vert:
            vertI_file = os.path.join(self.basePath + '_vertI.npy')
            vertIErr_file = os.path.join(self.basePath + '_vertIErr.npy')
            self.vertI = np.load(vertI_file)
            self.vertIErr = np.load(vertIErr_file)
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
            self.arrayParams[key] = value
        return True

    def get_array_param(self, key):
        return self.arrayParams[key]

    def get_running_param(self, key):
        return self.runningParams[key]

    def getNprime(self, VL, VR):
        Ch = np.array(self.get_array_param("Ch"))
        left_part_n_prime = np.copy(Ch[:, :-1])
        left_part_n_prime[:, 1:] = 0
        left_part_n_prime = flattenToRow(left_part_n_prime)
        right_part_n_prime = np.copy(Ch[:, 1:])
        right_part_n_prime[:, :-1] = 0
        right_part_n_prime = flattenToRow(right_part_n_prime)

        return self.n + (flattenToColumn(VL).dot(left_part_n_prime) + flattenToColumn(VR).dot(right_part_n_prime)).reshape(self.n.shape)

    def calc_hysteresis_score(self):
        IplusErr = self.I + self.IErr
        IminusErr = self.I - self.IErr
        score = -integrate.trapz(self.I, self.V)
        low_err = -integrate.trapz(IplusErr[:self.mid_idx], self.V[:self.mid_idx]) -\
                integrate.trapz(IminusErr[self.mid_idx:], self.V[self.mid_idx:])
        high_err = -integrate.trapz(IminusErr[:self.mid_idx], self.V[:self.mid_idx]) - \
                  integrate.trapz(IplusErr[self.mid_idx:], self.V[self.mid_idx:])
        return score, high_err-score, score-low_err

    def calc_blockade(self):
        I = self.I[:self.mid_idx]
        IErr = self.IErr[:self.mid_idx]
        matchingV = self.V[:self.mid_idx]
        score = np.min(np.abs(matchingV[I > 2*IErr]))
        return score, score, score

    def calc_jumps_score(self, window_size, up=True):
        score = 0
        high_err = 0
        low_err = 0
        if up:
            V = self.V[:self.mid_idx]
            I = self.I[:self.mid_idx]
            IplusErr = I + self.IErr[:self.mid_idx]
            IminusErr = I - self.IErr[:self.mid_idx]
        else:
            V = self.V[self.mid_idx:]
            I = self.I[self.mid_idx:]
            IplusErr = I + self.IErr[self.mid_idx:]
            IminusErr = I - self.IErr[self.mid_idx:]
        for v in V:
            in_window = np.logical_and(V >= v - window_size/2, V <= v + window_size/2)
            if in_window.any():
                new_score = np.max(I[in_window]) - np.min(I[in_window])
                if new_score > score:
                    score = new_score
                    high_err = np.max(IplusErr[in_window]) - np.min(IminusErr[in_window])
                    low_err = np.max(IminusErr[in_window]) - np.min(IplusErr[in_window])
        return score, high_err-score, score-low_err

    def calc_number_of_jumps(self):
        diffV1, diffV2, diff1, diff2, Vjumps1, Vjumps2 = self.calc_jumps_freq(self.I, self.IErr)
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

    def plot_jumps_freq(self, x_parameter, index=0, path=None):
        if x_parameter == "I":
            x = self.I
            xerr = self.IErr
        elif x_parameter == "n":
            n = self.getNprime(self.V, np.zeros(self.V.shape)).reshape(
                (self.n.shape[0], self.n.shape[1] * self.n.shape[2]))
            x = n[:,index]
            xerr = self.nErr[:,index]
            x/=500
            xerr/=500

        diffV1, diffV2, diff1, diff2, Vjumps1, Vjumps2 = self.calc_jumps_freq(x, xerr)
        plt.figure()
        plt.hist(np.hstack((diffV1,-diffV2)), bins=10)
        plt.xlabel('Delta V')
        plt.ylabel('Frequency')
        plt.title('Voltage between jumps histogram')
        if path is not None:
            plt.savefig(path + '_deltaV.png')
        plt.figure()
        plt.hist(np.hstack((diff1[diff1>0], -diff2[diff2<0])), bins=10)
        plt.xlabel('Delta I')
        plt.ylabel('Frequency')
        plt.title('Jumps heights histogram')
        if path is not None:
            plt.savefig(path + '_deltaI.png')
        plt.figure()
        plt.plot(self.V[:self.mid_idx], x[:self.mid_idx])
        for v in Vjumps1:
            marker = 'r--' if x_parameter == "I" else 'g--'
            plt.plot([v, v], [0, np.max(x)], marker)
        plt.plot(self.V[self.mid_idx:], x[self.mid_idx:])
        for v in Vjumps2:
            marker = 'm--' if x_parameter == "I" else 'c--'
            plt.plot([v, v], [0, np.max(x)], marker)
        plt.xlabel('V')
        plt.ylabel('I/n')
        # plt.xlim(2.35,2.5)
        if path is not None:
            plt.savefig(path + '_detected_jumps.png')

    def calc_jumps_freq(self, x, xerr):
        diff1 = np.diff(x[:self.mid_idx])
        diff2 = np.diff(x[self.mid_idx:])
        eps1 = np.clip(xerr[:self.mid_idx-1], a_min=0.001, a_max=None)
        eps2 = np.clip(xerr[self.mid_idx:-1], a_min=0.001, a_max=None)

        diff1[np.abs(diff1) < eps1] = 0
        diff2[np.abs(diff2) < eps2] = 0
        Vjumps1 = self.V[1:self.mid_idx]
        Vjumps1 = Vjumps1[diff1 > 0]
        Vjumps2 = self.V[self.mid_idx:-1]
        Vjumps2 = Vjumps2[diff2 < 0]
        xjumps1 = diff1[diff1 > 0]
        xjumps2 = diff2[diff2 < 0]
        new_Vjumps1 = []
        new_diff1 = []
        lastV = 0
        last_jump = 0
        for v,i in zip(Vjumps1, xjumps1):
            if v - lastV > 0.03:
                new_Vjumps1.append(v)
                new_diff1.append(last_jump + i)
                last_jump = 0
            else:
                last_jump += i
            lastV = v
        Vjumps1 = np.array(new_Vjumps1)
        diff1 = np.array(new_diff1)
        new_Vjumps2 = []
        new_diff2 = []
        lastV = np.max(self.V)
        last_jump = 0
        for v, i in zip(Vjumps2, xjumps2):
            if lastV - v > 0.03:
                new_Vjumps2.append(v)
                new_diff2.append(last_jump + i)
                last_jump = 0
            else:
                last_jump += i
            lastV = v
        Vjumps2 = np.array(new_Vjumps2)
        diff2 = np.array(new_diff2)
        diffV1 = np.diff(Vjumps1)
        diffV2 = np.diff(Vjumps2)
        return diffV1, diffV2, diff1, diff2, Vjumps1, Vjumps2


    def clac_fourier(self, eps=0.001):
        diff1 = np.diff(self.I[:self.mid_idx])
        diff2 = np.diff(self.I[self.mid_idx:])
        diff1[np.abs(diff1) < eps] = 0
        diff2[np.abs(diff2) < eps] = 0
        sample_space = self.V[1] - self.V[0]
        fou1 = np.abs(np.fft.rfft(diff1))
        fou2 = np.abs(np.fft.rfft(diff2))
        fou1[0]=0
        fou2[0]=0
        freq1 = np.fft.rfftfreq(diff1.size, sample_space)
        freq2 = np.fft.rfftfreq(diff2.size, sample_space)
        plt.figure()
        plt.plot(self.V[1:self.mid_idx], diff1, self.V[self.mid_idx:-1], diff2)
        plt.figure()
        plt.plot(freq1, fou1, label="increasing voltage")
        plt.plot(freq2, fou2, label="decreasing voltage")
        plt.xlabel("Frequency [1/V]")
        plt.ylabel("Fourier amplitude")
        plt.legend()

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


    def plot_array_params(self, parameter):
        M = self.runningParams["M"]*2 - 1
        N = self.runningParams["N"] + 1
        im = np.zeros(((M//2)*3+1,N*3))
        if parameter in ["R", "C"]:
            horzRows = np.arange(0, (M // 2) * 3 + 1, 3)
            horzCols = np.repeat(np.arange(0, 3 * N, 3), 2)
            horzCols[1::2] += 1
            vertRows = np.repeat(np.arange(1, (M // 2) * 3 + 1, 3), 2)
            vertRows[1::2] += 1
            vertCols = np.arange(2, 3 * N, 3)
            vertCols = vertCols[:-1]
            data_horz = np.array(self.arrayParams["Ch"]) if parameter == "C" else np.array(self.arrayParams["Rh"])
            data_vert = np.array(self.arrayParams["Cv"]) if parameter == "C" else np.array(self.arrayParams["Rv"])

            if parameter == "C":
                data_vert = data_vert[1:,:]
                cmap = 'Greens'
                neighbors = [data_horz[:, 1:], data_horz[:, :-1],
                             np.pad(data_vert, ((1, 0),(0,0))), np.pad(data_vert, ((0, 1),(0,0)))]
                data = np.average(neighbors, axis=0)
            elif parameter == "R":
                neighbors = [data_horz[:, 1:], data_horz[:, :-1],
                             np.pad(data_vert, ((1, 0),(0,0))), np.pad(data_vert, ((0, 1),(0,0)))]
                data = np.std(neighbors,axis=0)
                cmap = 'Reds'
            im[np.ix_(horzRows, horzCols)] = np.repeat(data_horz,2,axis=1)
            im[np.ix_(vertRows, vertCols)] = np.repeat(data_vert,2,axis=0)

        else:
            data = self.arrayParams["CG"] if parameter == "CG" else\
                self.arrayParams["RG"] if parameter == "RG" else self.arrayParams["VG"]
            cmap = 'RdBu' if parameter == "VG" else 'Greens'
        rows = np.arange(0, (M // 2) * 3 + 1, 3)
        cols = np.arange(2, 3 * N, 3)
        cols = cols[:-1]
        data = np.reshape(data, (self.runningParams["M"], self.runningParams["N"]))

        im[np.ix_(rows, cols)] = data
        data_max = np.max(im)
        data_min = np.min(im)
        plt.figure()
        plt.imshow(im, vmin=data_min, vmax=data_max, cmap=cmap,
                   aspect='equal')
        cb = plt.colorbar()
        cb.set_label(parameter)

    def createCapacitanceMatrix(self):
        """
        Creates the inverse capacitance matrix
        """
        Ch = np.array(self.get_array_param("Ch"))
        Cv = np.array(self.get_array_param("Cv"))
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
        n = self.getNprime(self.V, np.zeros(self.V.shape)).reshape((self.n.shape[0], self.n.shape[1] * self.n.shape[2]))
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
        ax.plot_surface(self.X, self.Y, Z)
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
        # if self.full:
        #     n = self.getNprime(self.V, np.zeros(self.V.shape)).reshape((self.n.shape[0], self.n.shape[1] * self.n.shape[2]))
        #     # Q = self.Q.reshape((self.Q.shape[0], self.Q.shape[1] * self.Q.shape[2]))
        #     # q = n + Q
        #     n_up = n[:self.mid_idx]
        #     n_down = n[self.mid_idx:]
        #     plt.figure()
        #     for i in range(20,30):
        #         plt.plot(V, np.flip(n_down[:,i]) - n_up[:,i])
        #     plt.ylabel("q")
        #     plt.xlabel("V")
        #     plt.xlim(xmin, np.max(V))



    def plot_results(self):
        plt.figure()
        plt.errorbar(self.V[:self.mid_idx], self.I[:self.mid_idx], fmt='r.', yerr=self.IErr[:self.mid_idx],
                     errorevery=10, label=self.directory + "_up")
        plt.errorbar(self.V[self.mid_idx:], self.I[self.mid_idx:], fmt='b.', yerr=self.IErr[self.mid_idx:],
                     errorevery=10, label=self.directory + "_down")
        plt.errorbar(self.alternativeV[:self.alternative_mid_idx], self.alternativeI[:self.alternative_mid_idx],
                     fmt='m.', yerr=self.alternativeIErr[:self.alternative_mid_idx],
                     errorevery=10, label=self.directory + "_alterntive")
        plt.errorbar(self.alternativeV[self.alternative_mid_idx:], self.alternativeI[self.alternative_mid_idx:],
                     fmt='c.', yerr=self.alternativeIErr[self.alternative_mid_idx:],
                     errorevery=10, label=self.directory + "_alterntive")
        plt.xlabel('Voltage')
        plt.ylabel('Current')
        plt.legend()
        if self.full:
            n = self.getNprime(self.V, np.zeros(self.V.shape)).reshape((self.n.shape[0], self.n.shape[1] * self.n.shape[2]))
            Q = self.Q.reshape((self.Q.shape[0], self.Q.shape[1] * self.Q.shape[2]))
            q = n+Q
            nErr = self.nErr.reshape((self.nErr.shape[0], self.nErr.shape[1] *self. nErr.shape[2]))
            QErr = self.QErr.reshape((self.QErr.shape[0], self.QErr.shape[1] * self.QErr.shape[2]))
            nplusErr = n + nErr
            nminusErr = n - nErr
            QplusErr = Q + QErr
            QminusErr = Q - QErr
            plt.figure()
            for i in range(len(n[0])):
                plt.plot(self.V[:self.mid_idx], n[:self.mid_idx, i], 'b',
                         self.V[self.mid_idx:], n[self.mid_idx:, i], 'r')
                plt.xlabel('Voltage')
                plt.ylabel('Occupation')
            # factor = np.max(np.diff(np.sum(n[:self.mid_idx, :],axis=1)))/np.max(np.diff(self.I[:self.mid_idx]),)
            # Idiff = np.diff(self.I[:self.mid_idx])
            # Idiff[Idiff < 0.0001] = 0
            # plt.plot(self.V[1:self.mid_idx], np.diff(np.sum(n[:self.mid_idx, :],axis=1)), 'g',
            #          self.V[1:self.mid_idx], Idiff, 'orange')
            # plt.xlabel('Voltage')
            # plt.ylabel('Occupation')
            # # plt.figure()
            # plt.plot(self.V[:self.mid_idx], np.sum(n[:self.mid_idx, :],axis=1)/40000, 'b',
            #          self.V[self.mid_idx:], np.sum(n[self.mid_idx:, :],axis=1)/40000, 'r')
            # plt.xlabel('Voltage')
            # plt.ylabel('Total Occupation')
            plt.figure()
            for i in range(len(Q[0])):
                plt.plot(self.V[:self.mid_idx], -Q[:self.mid_idx, i], 'g',
                         self.V[self.mid_idx:], -Q[self.mid_idx:, i], 'c',
                         self.V[:self.mid_idx], -QplusErr[:self.mid_idx, i], 'g--',
                         self.V[self.mid_idx:], -QplusErr[self.mid_idx:, i], 'c--',
                         self.V[:self.mid_idx], -QminusErr[:self.mid_idx, i], 'g--',
                         self.V[self.mid_idx:], -QminusErr[self.mid_idx:, i], 'c--')
                plt.xlabel('Voltage')
                plt.ylabel('Chagre')

            plt.figure()
            for i in range(len(self.full_I)):
                plt.plot(self.V[:self.mid_idx], self.full_I[i,:self.mid_idx], 'o',
                        self.V[self.mid_idx:], self.full_I[i,self.mid_idx:], '*')
                plt.xlabel('Voltage')
                plt.ylabel('Current')

class MultiResultAnalyzer:
    """ Used for statistical analysis of results from many simulations"""
    def __init__(self, directories_list, files_list, relevant_running_params=None, relevant_array_params=None, out_directory = None,
                 groups=None, group_names=None, resistance_line=0, full=False):
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
        self.hysteresisScores = []
        self.hysteresisScoresHighErr = []
        self.hysteresisScoresLowErr = []
        self.jumpScores = []
        self.jumpScoresHighErr = []
        self.jumpScoresLowErr = []
        self.blockadeScores = []
        self.blockadeScoresHighErr = []
        self.blockadeScoresLowErr = []
        self.jumpsNum = []
        self.resistanceScore = []
        if relevant_array_params:
            self.arrayParams = {p: [] for p in relevant_array_params}
        else:
            self.arrayParams = dict()
        if relevant_running_params:
            self.runningParams = {p: [] for p in relevant_running_params}
        else:
            self.runningParams = dict()
        self.disorders = {p: [] for p in ["C","R","CG","VG"]}
        self.groups = np.array(groups)
        self.groupNames = group_names
        self.load_data(resistance_line, full=full)

    def load_data(self, line=0, full=False):
        """
        Loading results and calculating scores
        """
        for directory, fileName in zip(self.directories, self.fileNames):
            processor = SingleResultsProcessor(directory, fileName, fullOutput=full)
            processor.load_results()
            hyst, hystHigh, hystLow = processor.calc_hysteresis_score()
            block, blockHigh, blockLow = processor.calc_blockade()
            jump, jumpHigh, jumpLow = processor.calc_jumps_score(0.1)
            num_jumps = processor.calc_number_of_jumps()
            resistance_score = processor.calc_resistance_hysteresis_score(line)
            self.resistanceScore.append(resistance_score)
            self.jumpsNum.append(num_jumps)
            self.hysteresisScores.append(hyst)
            self.hysteresisScoresHighErr.append(hystHigh)
            self.hysteresisScoresLowErr.append(hystLow)
            self.blockadeScores.append(block)
            self.blockadeScoresHighErr.append(blockHigh)
            self.blockadeScoresLowErr.append(blockLow)
            self.jumpScores.append(jump)
            self.jumpScoresHighErr.append(jumpHigh)
            self.jumpScoresLowErr.append(jumpLow)
            processor.load_params()
            for p in self.arrayParams:
                self.arrayParams[p].append(processor.get_array_param(p))
            for p in self.runningParams:
                self.runningParams[p].append(processor.get_running_param(p))
            for p in self.disorders:
                self.disorders[p].append(self.calc_disorder(p, processor))

    def calc_disorder(self, name, processor):
        if name in ["C","R"]:
            var = np.var(processor.get_array_param(name + "h")) + np.var(processor.get_array_param(name + "v"))
            std = np.sqrt(var)
        else:
            std = np.std(processor.get_array_param(name))
        return std

    def get_relevant_scores(self, scoreName):
        if scoreName == 'hysteresis':
            scores = np.array(self.hysteresisScores)
            scoresHighErr = np.array(self.hysteresisScoresHighErr)
            scoresLowErr = np.array(self.hysteresisScoresLowErr)
        elif scoreName == 'jump':
            scores = np.array(self.jumpScores)
            scoresHighErr = np.array(self.jumpScoresHighErr)
            scoresLowErr = np.array(self.jumpScoresLowErr)
        elif scoreName == 'blockade':
            scores = np.array(self.blockadeScores)
            scoresHighErr = np.array(self.blockadeScoresHighErr)
            scoresLowErr = np.array(self.blockadeScoresLowErr)
        elif scoreName == 'jumpsNum':
            scores = np.array(self.jumpsNum)
            scoresHighErr = scores
            scoresLowErr = scores
        elif scoreName == 'resistance':
            scores = np.array(self.resistanceScore)
            scoresHighErr = scores
            scoresLowErr = scores
        else:
            raise NotImplementedError
        return scores, scoresHighErr, scoresLowErr


    def plot_hystogram(self, scoreName, runinng_parameter_vals, array_parameter_vals, title):
        """
        Plotting hystogram of the given score for results with parameters matching the given values
        :param scoreName: name of score (hysteresis, jump, blockade)
        :param (running/array)_parameter_vals: dictionary with {parameter_name:list_of_parameter_values}.
        For non-given parameters all values will be considered as valid
        """
        scores, _, _ = self.get_relevant_scores(scoreName)
        indices = set(range(scores.size))
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
        scores = scores[indices]
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
        score1, _, _ = self.get_relevant_scores(score1_name)
        score2, _, _ = self.get_relevant_scores(score2_name)
        plt.figure()
        for i in range(np.max(self.groups)+1):
            plt.scatter(score1[self.groups == i],score2[self.groups == i], label=self.groupNames[i])
        plt.legend()
        plt.xlabel(score1_name)
        plt.ylabel(score2_name)
        plt.savefig(os.path.join(self.outDir, score1_name + "_vs_" + score2_name + '.png'))



    def get_y_label(self, score_name):
        if score_name == "hysteresis":
            return "loop area"
        elif score_name == "jump":
            return "voltage diff for biggest jump"
        elif score_name == "blockade":
            return "threshold voltage"


    def plot_score(self, score_name):
        plt.figure(figsize=(16,18))
        y, y_high_err, y_low_err = self.get_relevant_scores(score_name)
        plt.plot(np.arange(y.size), y)
        plt.ylabel(score_name + " score")
        plt.xlabel("example num.")


    def plot_score_by_parameter(self, score_name, parameter_names, title, runningParam=True, plot=True):
        """
        Plots score as a function of given parameter
        :param score_name: relevant score name (hysteresis, jump, blockade)
        :param parameter_name: relevant parameter
        """
        if len(parameter_names) == 1:
            parameter_name = parameter_names[0]
            y, y_high_err, y_low_err = self.get_relevant_scores(score_name)
            if runningParam:
                x = self.runningParams[parameter_name]
            else:
                x = self.arrayParams[parameter_name]
            x, y, y_high_err, y_low_err = self.average_similar_results(x, y, y_high_err, y_low_err)
            errors = [y_low_err, y_high_err]
            fig = plt.figure()
            plt.errorbar(x, y, yerr=errors, marker='o')
            plt.title(title)
            plt.xlabel(parameter_name)
            plt.ylabel(self.get_y_label(score_name))
            plt.savefig(os.path.join(self.outDir, title.replace(' ', '_') + '.png'))
            plt.close(fig)
            np.save(os.path.join(self.outDir, title + "parameter"), np.array(x))
            np.save(os.path.join(self.outDir, title + "score"), np.array(y))
            np.save(os.path.join(self.outDir, title + "score_high_err"), np.array(y_high_err))
            np.save(os.path.join(self.outDir, title + "score_low_err"), np.array(y_low_err))
        elif len(parameter_names) == 2:
            z, z_high_err, z_low_err = self.get_relevant_scores(score_name)
            if runningParam:
                x = self.runningParams[parameter_names[0]]
                y = self.runningParams[parameter_names[1]]
            else:
                x = self.arrayParams[parameter_names[0]]
                y = self.arrayParams[parameter_names[1]]
            points = [(x[i],y[i]) for i in range(len(x))]
            points, z, z_high_err, z_low_err = self.average_similar_results(points, z, z_high_err, z_low_err)
            x = [point[0] for point in points]
            y = [point[1] for point in points]
            fig = plt.figure()
            plt.scatter(x, y, c=z, cmap='viridis',norm=colors.LogNorm())
            plt.title(title)
            plt.xlabel(parameter_names[0])
            plt.ylabel(parameter_names[1])
            cbar = plt.colorbar()
            cbar.set_label(self.get_y_label(score_name))
            plt.savefig(os.path.join(self.outDir, title.replace(' ', '_') + '.png'))
            plt.close(fig)
        else:
            raise NotImplementedError

    def plot_results_by_disorder(self, parameters, scores):
        for param in parameters:
            for score in scores:
                x = self.disorders[param]
                y, y_high_err, y_low_err = self.get_relevant_scores(score)
                plt.figure()
                plt.scatter(x,y)
                plt.xlabel(param + " disorder")
                plt.ylabel(score)
                if self.outDir is not None:
                    plt.savefig(os.path.join(self.outDir, score+"_"+ param +'_'+'disorder.png'))

    def average_similar_results(self, xs, ys, ys_high_err, ys_low_err):
        new_x = []
        new_y = []
        new_y_high_err = []
        new_y_low_err = []
        for x in xs:
            if x not in new_x:
                new_x.append(x)
                indices = [i for i, val in enumerate(xs) if val == x]
                relevant_ys = np.array(ys[indices])
                relevant_high_err = np.array(ys_high_err[indices])
                relevant_low_err = np.array(ys_low_err[indices])
                new_y.append(np.average(relevant_ys))
                new_y_high_err.append(np.sqrt(np.average(relevant_high_err**2)/relevant_high_err.size))
                new_y_low_err.append(np.sqrt(np.average(relevant_low_err**2)/relevant_low_err.size))
        return new_x, new_y, new_y_high_err, new_y_low_err

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



def bistabilityAnalysis(x):
    avg = np.average(x)
    var = np.var(x)
    dist_from_mean = (x - avg)**2
    bins = [x]
    if np.max(dist_from_mean) > var:
        model = GaussianMixture(n_components=2)
        labels = model.fit_predict(x.reshape(-1,1))
        bins = [x[labels==0], x[labels==1]]
    return bins




if __name__ == "__main__":
    if len(sys.argv) > 1:
        action = sys.argv[1]
    else:
        print("No action was selected")
        exit(0)
    if action == "perpendicular":
        ###### Perpendicular Current analysis ######
        directory = "2d_long_array_bgu_with_perp_point_contact"
        names = ["array_5_15_r_std_9_run_" + str(run) for run in range(1,11)]
        full = True
        save_re_analysis = True
        for name in names:
            s = SingleResultsProcessor(directory, name, fullOutput=full, vertCurrent=True)
            if save_re_analysis:
                s.save_re_analysis()
            s.plot_resistance(out=os.path.join(directory,name))

    elif action == "jumps":
        ###### Jumps analysis ############
        directory = "bgu_2d_narrow_arrays"
        names = ["array_5_15_disorder_c_std_0.1_r_std_9_run_10", "array_5_15_disorder_c_std_0_r_std_9_run_10",
                 "array_5_15_disorder_c_std_0.1_r_std_9_run_10_high_r_disorder"]
        full = True
        save_re_analysis = False
        for name in names:
            s = SingleResultsProcessor(directory, name, fullOutput=full, vertCurrent=False)
            if save_re_analysis:
                s.save_re_analysis()
            s.plot_jumps_freq("I")
            # s.plot_results()
        plt.show()
    else:
        ###### General analysis ############
        directory = "bgu_2d_narrow_arrays"
        names = ["array_5_15_disorder_c_std_0.1_r_std_9_run_4", "array_5_15_disorder_c_std_1_r_std_9_run_4",
                 "array_5_15_disorder_c_std_1.5_r_std_9_run_4"]
        full = True
        save_re_analysis = False
        for name in names:
            s = SingleResultsProcessor(directory, name, fullOutput=full, vertCurrent=False)
            if save_re_analysis:
                s.save_re_analysis()
            # s.plot_jumps_freq("I")
            s.plot_results()
            s.plot_array_params("R")
            s.plot_array_params("C")
        plt.show()
