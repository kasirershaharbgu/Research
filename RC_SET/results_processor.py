import numpy as np
from matplotlib import pyplot as plt
from scipy import integrate
from ast import literal_eval
import os
import re

EPS = 1e-6
class SingleResultsProcessor:
    """ Used for proccessing result from one dot array simulation"""
    def __init__(self, directory, file_name, fullOutput=False):
        self.basePath = os.path.join(directory, file_name)
        self.fileName = file_name
        self.directory = directory
        self.full = fullOutput
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
        self.load_results(fullOutput=fullOutput)
        self.load_params()

    def load_results(self, fullOutput=False):
        I_file = os.path.join(self.basePath + '_I.npy')
        V_file = os.path.join(self.basePath + '_V.npy')
        IErr_file = os.path.join(self.basePath + '_IErr.npy')
        self.I = np.load(I_file)
        self.IErr = np.load(IErr_file)
        # self.IErr = np.zeros(self.I.shape)
        self.V = np.load(V_file)
        self.mid_idx = self.V.size // 2
        if fullOutput:
            n_file = os.path.join(self.basePath + '_n.npy')
            Q_file = os.path.join(self.basePath + '_Q.npy')
            nErr_file = os.path.join(self.basePath + '_nErr.npy')
            QErr_file = os.path.join(self.basePath + '_QErr.npy')
            self.n = np.load(n_file)
            self.Q = np.load(Q_file)
            self.nErr = np.load(nErr_file)
            self.QErr = np.load(QErr_file)
            # self.nErr = np.zeros(self.n.shape)
            # self.QErr = np.zeros(self.Q.shape)
        return True

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
        IminusErr = self.I[:self.mid_idx] - self.IErr[:self.mid_idx]
        IplusErr = self.I[:self.mid_idx] + self.IErr[:self.mid_idx]
        matchingV = self.V[:self.mid_idx]
        score = np.min(np.abs(matchingV[I > EPS]))
        high_err = np.min(np.abs(matchingV[IminusErr > EPS]))
        low_err = np.min(np.abs(matchingV[IplusErr > EPS]))
        return score, high_err-score, score-low_err

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

    def plot_results(self):
        IplusErr = self.I + self.IErr
        IminusErr = self.I - self.IErr
        plt.figure()
        plt.plot(self.V[:self.mid_idx], self.I[:self.mid_idx], 'b',
                 self.V[self.mid_idx:], self.I[self.mid_idx:], 'r',
                 self.V[:self.mid_idx], IplusErr[:self.mid_idx], 'b--',
                 self.V[self.mid_idx:], IplusErr[self.mid_idx:],'r--',
                 self.V[:self.mid_idx], IminusErr[:self.mid_idx], 'b--',
                 self.V[self.mid_idx:], IminusErr[self.mid_idx:], 'r--')
        plt.xlabel('Voltage')
        plt.ylabel('Current')
        if self.full:
            plt.figure()
            n = self.n.reshape((self.n.shape[0], self.n.shape[1] * self.n.shape[2]))
            Q = self.Q.reshape((self.Q.shape[0], self.Q.shape[1] * self.Q.shape[2]))
            nErr = self.nErr.reshape((self.nErr.shape[0], self.nErr.shape[1] *self. nErr.shape[2]))
            QErr = self.QErr.reshape((self.QErr.shape[0], self.QErr.shape[1] * self.QErr.shape[2]))
            nplusErr = n + nErr
            nminusErr = n - nErr
            QplusErr = Q + QErr
            QminusErr = Q - QErr
            for i in range(len(n[0])):
                plt.plot(self.V[:self.mid_idx], n[:self.mid_idx, i], 'b',
                         self.V[self.mid_idx:], n[self.mid_idx:, i], 'r',
                         self.V[:self.mid_idx], nplusErr[:self.mid_idx, i], 'b--',
                         self.V[self.mid_idx:], nplusErr[self.mid_idx:, i], 'r--',
                         self.V[:self.mid_idx], nminusErr[:self.mid_idx, i], 'b--',
                         self.V[self.mid_idx:], nminusErr[self.mid_idx:, i], 'r--')
                plt.xlabel('Voltage')
                plt.ylabel('Occupation')
            plt.figure()
            for i in range(len(Q[0])):
                plt.plot(self.V[:self.mid_idx], Q[:self.mid_idx, i], 'b',
                         self.V[self.mid_idx:], Q[self.mid_idx:, i], 'r',
                         self.V[:self.mid_idx], QplusErr[:self.mid_idx, i], 'b--',
                         self.V[self.mid_idx:], QplusErr[self.mid_idx:, i], 'r--',
                         self.V[:self.mid_idx], QminusErr[:self.mid_idx, i], 'b--',
                         self.V[self.mid_idx:], QminusErr[self.mid_idx:, i], 'r--')
                plt.xlabel('Voltage')
                plt.ylabel('Chagre')

class MultiResultAnalyzer:
    """ Used for statistical analysis of results from many simulations"""
    def __init__(self, directories_list, files_list, relevant_running_params, relevant_array_params, out_directory):
        """
        Initializing analyzer
        :param directories_list: list of directories for result files
        :param files_list: list of result names
        :param relevant_array_params: list of array parameters that are relevant
        :param relevant_running_params: list of array parameters that are relevant
        :param out_directory: directory for results
        """
        self.outDir = out_directory
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
        self.arrayParams = {p: [] for p in relevant_array_params}
        self.runningParams = {p: [] for p in relevant_running_params}
        self.load_data()

    def load_data(self):
        """
        Loading results and calculating scores
        """
        for directory, fileName in zip(self.directories, self.fileNames):
            processor = SingleResultsProcessor(directory, fileName)
            processor.load_results()
            hyst, hystHigh, hystLow = processor.calc_hysteresis_score()
            block, blockHigh, blockLow = processor.calc_blockade()
            jump, jumpHigh, jumpLow = processor.calc_jumps_score(0.1)
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
        plt.savefig(os.path.join(self.outDir, title.replace(' ', '_') + '.png'))
        plt.close(fig)
        np.save(os.path.join(self.outDir,title), np.array(scores))

    def get_y_label(self, score_name):
        if score_name == "hysteresis":
            return "loop area"
        elif score_name == "jump":
            return "voltage diff for biggest jump"
        elif score_name == "blockade":
            return "threshold voltage"

    def plot_score(self, score_name, parameter_name, title, runningParam=True):
        """
        Plots score asa function of given parameter
        :param score_name: relevant score name (hysteresis, jump, blockade)
        :param parameter_name: relevant parameter
        """
        y, y_high_err, y_low_err = self.get_relevant_scores(score_name)
        if runningParam:
            x = self.runningParams[parameter_name]
        else:
            x = self.runningParams[parameter_name]
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

    def average_similar_results(self, xs, ys, ys_high_err, ys_low_err):
        new_x = []
        new_y = []
        new_y_high_err = []
        new_y_low_err = []
        for x in xs:
            if x not in new_x:
                new_x.append(x)
                indices = np.array(xs) == x
                relevant_ys = np.array(ys[indices])
                relevant_high_err = np.array(ys_high_err[indices])
                relevant_low_err = np.array(ys_low_err[indices])
                new_y.append(np.average(relevant_ys))
                new_y_high_err.append(np.sqrt(np.average(relevant_high_err**2)/relevant_high_err.size))
                new_y_low_err.append(np.sqrt(np.average(relevant_low_err**2)/relevant_low_err.size))
        return new_x, new_y, new_y_high_err, new_y_low_err


if __name__ == "__main__":
    # directory_list = ["C:\\Users\\shahar\\Research\\RC_SET\\3X3_array_statistics_r_avg_10"] * 90
    # files_list = ["c_std_0.5_r_std_" + str(i) + "_run_" + str(j) for i in range(1,10) for j in range(1,11)]
    # m = MultiResultAnalyzer(directory_list, files_list, ["C_std", "R_std"], [],
    #                         "C:\\Users\\shahar\\Research\\old_results\\3X3_array_statistics_r_avg_10")
    # m.plot_score('hysteresis', 'R_std', 'c_std_0.5_hysteresis')
    # m.plot_score('jump', 'R_std', 'c_std_0.5_jump')
    # m.plot_score('blockade', 'R_std', 'c_std_0.5_blockade')
    # for r_std in range(1,10):
    #     for score in ['hysteresis', 'jump', 'blockade']:
    #         m.plot_hystogram(score, {"R_std": [r_std], "C_std": [0.5]}, {},
    #                          "c_std_0.5_r_std_" + str(r_std) + "_" + score + "_hystogram")
    # for score in ['hysteresis', 'jump', 'blockade']:
    #     m.plot_hystogram(score, {"C_std": [0.4]}, {},
    #                          "c_std_0.5_" + score + "_hystogram")
    s = SingleResultsProcessor("finit_dos_shorter_run","array_3_3_finit_dos_cg_disorder_refine",fullOutput=True)
    s.plot_results()
    plt.show()





