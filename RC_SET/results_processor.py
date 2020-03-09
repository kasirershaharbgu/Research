import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors
from scipy import integrate
from ast import literal_eval
import os
import re
from sklearn.mixture import GaussianMixture
from mpl_toolkits.mplot3d import Axes3D

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

    def reAnalyzeI(self):
        full_I = self.full_I.T
        for i,point in enumerate(full_I):
            bins = bistabilityAnalysis(point)
            if i > 0 and len(bins) > 1 and len(bins[0]) > 0 and len(bins[1]) > 0:
                avg1 = np.average(bins[0])
                avg2 = np.average(bins[1])
                if np.abs(avg1 - self.I[i-1]) < np.abs(avg2 - self.I[i-1]):
                # if len(bins[0]) > len(bins[1]):
                    self.I[i] = avg1
                    self.IErr[i] = np.std(bins[0])
                else:
                    self.I[i] = avg2
                    self.IErr[i] = np.std(bins[1])
            else:
                if len(bins) > 1:
                    bin = bins[0] if len(bins[0]) > len(bins[1]) else bins[1]
                else:
                    bin = bins[0]
                self.I[i] = np.average(bin)
                self.IErr[i] = np.std(bin)
        return True


    def load_results(self, fullOutput=False):
        I_file = os.path.join(self.basePath + '_I.npy')
        V_file = os.path.join(self.basePath + '_V.npy')
        IErr_file = os.path.join(self.basePath + '_IErr.npy')
        self.I = np.load(I_file)
        self.IErr = np.load(IErr_file)
        self.V = np.load(V_file)
        self.mid_idx = self.V.size // 2
        if fullOutput:
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
        return True

    def save_re_analysis(self):
        np.save(self.basePath + '_I_clean', self.I)
        np.save(self.basePath + '_I_Err_clean', self.IErr)
        IplusErr = self.I + self.IErr
        IminusErr = self.I - self.IErr
        fig = plt.figure()
        plt.plot(self.V[:self.mid_idx], self.I[:self.mid_idx], 'b.',
                 self.V[self.mid_idx:], self.I[self.mid_idx:], 'r.',
                 self.V[:self.mid_idx], IplusErr[:self.mid_idx], 'b--',
                 self.V[self.mid_idx:], IplusErr[self.mid_idx:], 'r--',
                 self.V[:self.mid_idx], IminusErr[:self.mid_idx], 'b--',
                 self.V[self.mid_idx:], IminusErr[self.mid_idx:], 'r--')
        plt.xlabel('Voltage')
        plt.ylabel('Current')
        plt.savefig(self.basePath + '_IV_clean')
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

    def calc_jumps_freq(self):
        diff1 = np.diff(self.I[:self.mid_idx])
        diff2 = np.diff(self.I[self.mid_idx:])
        diff1[np.abs(diff1) < 0.002] = 0
        diff2[np.abs(diff2) < 0.002] = 0
        Vjumps1 = self.V[1:self.mid_idx]
        Vjumps1 = Vjumps1[diff1 > 0]
        Vjumps2 = self.V[self.mid_idx:-1]
        Vjumps2 = Vjumps2[diff2 < 0]
        diffV1 = np.diff(Vjumps1)
        diffV2 = np.diff(Vjumps2)
        plt.figure()
        plt.hist(diffV1, bins=10)
        plt.figure()
        plt.hist(-diffV2, bins=10)
        plt.figure()
        plt.hist(diff1[diff1>0], bins=10)
        plt.figure()
        plt.hist(-diff2[diff2<0], bins=10)

    def clac_fourier(self):
        diff1 = np.diff(self.I[:self.mid_idx])
        diff2 = np.diff(self.I[self.mid_idx:])
        diff1[np.abs(diff1) < 0.002] = 0
        diff2[np.abs(diff2) < 0.002] = 0
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
        plt.plot(freq1, fou1, freq2, fou2)


    def plot_results(self):
        # IplusErr = self.I + self.IErr
        # IminusErr = self.I - self.IErr
        # plt.figure()
        # plt.plot(self.V[:self.mid_idx], self.I[:self.mid_idx], 'b.',
        #          self.V[self.mid_idx:], self.I[self.mid_idx:], 'r.',
        #          self.V[:self.mid_idx], IplusErr[:self.mid_idx], 'b--',
        #          self.V[self.mid_idx:], IplusErr[self.mid_idx:],'r--',
        #          self.V[:self.mid_idx], IminusErr[:self.mid_idx], 'b--',
        #          self.V[self.mid_idx:], IminusErr[self.mid_idx:], 'r--')
        # plt.xlabel('Voltage')
        # plt.ylabel('Current')
        if self.full:
            # plt.figure()
            n = self.n.reshape((self.n.shape[0], self.n.shape[1] * self.n.shape[2]))
            Q = self.Q.reshape((self.Q.shape[0], self.Q.shape[1] * self.Q.shape[2]))
            q = n+Q
            nErr = self.nErr.reshape((self.nErr.shape[0], self.nErr.shape[1] *self. nErr.shape[2]))
            QErr = self.QErr.reshape((self.QErr.shape[0], self.QErr.shape[1] * self.QErr.shape[2]))
            nplusErr = n + nErr
            nminusErr = n - nErr
            QplusErr = Q + QErr
            QminusErr = Q - QErr
            for i in range(6,100,10):
                plt.figure()
                plt.plot(self.V[:self.mid_idx], n[:self.mid_idx, i], 'b',
                         self.V[self.mid_idx:], n[self.mid_idx:, i], 'r',
                         self.V[:self.mid_idx], nplusErr[:self.mid_idx, i], 'b--',
                         self.V[self.mid_idx:], nplusErr[self.mid_idx:, i], 'r--',
                         self.V[:self.mid_idx], nminusErr[:self.mid_idx, i], 'b--',
                         self.V[self.mid_idx:], nminusErr[self.mid_idx:, i], 'r--')
                plt.xlabel('Voltage')
                plt.ylabel('Occupation')
            # plt.figure()
            # plt.plot(self.V[:self.mid_idx], np.sum(n[:self.mid_idx, :],axis=1), 'b',
            #          self.V[self.mid_idx:], np.sum(n[self.mid_idx:, :],axis=1), 'r')
            # plt.xlabel('Voltage')
            # plt.ylabel('Total Occupation')
            # plt.figure()
            # for i in range(len(Q[0])):
            #     plt.plot(self.V[:self.mid_idx], Q[:self.mid_idx, i], 'b',
            #              self.V[self.mid_idx:], Q[self.mid_idx:, i], 'r',
            #              self.V[:self.mid_idx], QplusErr[:self.mid_idx, i], 'b--',
            #              self.V[self.mid_idx:], QplusErr[self.mid_idx:, i], 'r--',
            #              self.V[:self.mid_idx], QminusErr[:self.mid_idx, i], 'b--',
            #              self.V[self.mid_idx:], QminusErr[self.mid_idx:, i], 'r--')
            #     plt.xlabel('Voltage')
            #     plt.ylabel('Chagre')
            # plt.figure()
            # for i in range(len(Q[0])):
            #     plt.plot(self.V[:self.mid_idx], q[:self.mid_idx, i], 'b',
            #              self.V[self.mid_idx:], q[self.mid_idx:, i], 'r')
            #     plt.xlabel('Voltage')
            #     plt.ylabel('Chagre on tunneling junctions')
            # plt.figure()
            # for i in range(len(self.full_I)):
            #     plt.plot(self.V[:self.mid_idx], self.full_I[i,:self.mid_idx], '.',
            #             self.V[self.mid_idx:], self.full_I[i,self.mid_idx:], '*')
            #     plt.xlabel('Voltage')
            #     plt.ylabel('Chagre')

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

    def plot_score(self, score_name, parameter_names, title, runningParam=True, plot=True):
        """
        Plots score asa function of given parameter
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
    # for c_std in [0.1, 0.2, 0.3, 0.4, 0.5]:
    #     files_list = ["c_std_" + str(c_std) + "_r_std_" + str(i) + "_run_" + str(j) for i in range(1,10) for j in range(1,11)]
    #     directory_list = ["C:\\Users\\shahar\\Research\\old_results\\3X3_array_statistics_r_avg_10"] * len(files_list)
    #     m = MultiResultAnalyzer(directory_list, files_list, ["C_std", "R_std"], [],
    #                             "C:\\Users\\shahar\\Research\\old_results\\3X3_array_statistics_r_avg_10")
    #     m.plot_score('hysteresis', ['R_std'], 'c_std_' + str(c_std) + '_hysteresis')
    #     m.plot_score('jump', ['R_std'], 'c_std_' + str(c_std) + '_jump')
    #     m.plot_score('blockade', ['R_std'], 'c_std_' + str(c_std) + '_blockade')
    #     for r_std in range(1,10):
    #         for score in ['hysteresis', 'jump', 'blockade']:
    #             m.plot_hystogram(score, {"R_std": [r_std], "C_std": [c_std]}, {},
    #                              "c_std_" + str(c_std) + "_r_std_" + str(r_std) + "_" + score + "_hystogram")
    #     for score in ['hysteresis', 'jump', 'blockade']:
    #         m.plot_hystogram(score, {"C_std": [c_std]}, {},
    #                              "c_std_" +  str(c_std) + "_" + score + "_hystogram")

    # files_list = ["c_std_" + str(c_std) + "_r_std_" + str(i) + "_run_" + str(j) for c_std in [0.1,0.2,0.3,0.4,0.5] for i in range(1,10) for j in range(1,11)]
    # directory_list = ["C:\\Users\\shahar\\Research\\old_results\\3X3_array_statistics_r_avg_10"] * len(files_list)
    # m = MultiResultAnalyzer(directory_list, files_list, ["C_std", "R_std"], [],
    #                         "C:\\Users\\shahar\\Research\\old_results\\3X3_array_statistics_r_avg_10")
    # m.plot_score('hysteresis', ['R_std','C_std'], 'all_hysteresis')
    # m.plot_score('jump', ['R_std','C_std'], 'all_jump')
    # m.plot_score('blockade', ['R_std','C_std'], 'all_blockade')

    # directory = "2d_array_bgu"
    # for name in ['array_10_10_r_disorder_cg_disorder_run_3',
    #              'array_10_10_r_disorder_cg_disorder_variable_ef_run_2',
    #              'array_10_10_r_disorder_cg_disorder_variable_ef_run_3',
    #              'array_10_10_r_disorder_run_1',
    #              'array_10_10_r_disorder_run_3']:
    #     s = SingleResultsProcessor(directory, name,
    #                                    fullOutput=True)
    #     s.save_re_analysis()

    directory = "2d_array_bgu"
    name = "array_10_10_r_disorder_run_"
    for run in ["1"]:
        s = SingleResultsProcessor(directory, name+run,fullOutput=True)
        # s.calc_jumps_freq()
        # s.clac_fourier()
        s.plot_results()
    plt.show()





