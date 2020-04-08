__author__ = 'shahar'
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from scipy.linalg import solve as linsolve
from scipy.integrate import odeint
import os
from optparse import OptionParser

class ThermalArray:
    """
    A model for array with different local temperature for each location
    """
    def __init__(self, rows, columns, density_avg, density_std, resistance,
                 diffusionConst, bath_temp, VL, VR, VU, VD):
        self.rows = rows
        self.columns = columns
        self.diffusionConst = diffusionConst
        self.R0 = resistance
        self.VL = VL
        self.VR = VR
        self.VU = VU
        self.VD = VD
        self.tempConst = create_random_array(rows, columns, navg, nstd,
                                             "exp", True)
        self.bath_temperature = bath_temp
        self.setLaplacianMatrix()
    def getResistance(self, T):
        return self.R0*np.exp(1/T)

    def getConductanceMatrix(self, R):
        padded = np.pad(R, ((1,1),(1,1)), mode='constant')
        diagonal = 1/(R + padded[2:,1:-1]) +\
                   1/(R + padded[:-2,1:-1]) +\
                   1/(R + padded[1:-1, 2:]) +\
                   1/(R + padded[1:-1, :-2])
        second_diagonal = 1/(R + padded[1:-1,:-2])
        second_diagonal[:,0] = 0
        second_diagonal = second_diagonal.flatten()
        second_diagonal = second_diagonal[1:]
        n_diagonal = 1/(R[1:,:]  + R[:-1,:])
        return np.diagflat(diagonal) - \
               np.diagflat(second_diagonal, k=-1) - \
               np.diagflat(second_diagonal, k=1) - \
               np.diagflat(n_diagonal, k=-self.columns) - \
               np.diagflat(n_diagonal, k=self.columns)

    def setLaplacianMatrix(self):
        initial = np.ones((self.rows, self.columns))
        padded = np.pad(initial, ((1, 1), (1, 1)), mode='constant')
        diagonal = padded[2:, 1:-1] + padded[:-2, 1:-1] + \
                   padded[1:-1, 2:] + padded[1:-1, :-2]
        second_lower_diagonal = padded[1:-1, :-2].flatten()
        second_upper_diagonal = padded[1:-1, 2:].flatten()
        n_upper_diagonal = initial[1:, :]
        n_lower_diagonal = initial[:-1, :]
        self.laplacian = -np.diagflat(diagonal) + \
                     np.diagflat(second_lower_diagonal[1:], k=-1) + \
                     np.diagflat(second_upper_diagonal[:-1], k=1) + \
                     np.diagflat(n_lower_diagonal, k=-self.columns) + \
                     np.diagflat(n_upper_diagonal, k=self.columns)

    def getEdgeVec(self, R):
        edgeMask = np.zeros(R.shape)
        edgeMask[0,:] += self.VU
        edgeMask[-1,:] += self.VD
        edgeMask[:,0] += self.VL
        edgeMask[:,-1] += self.VR
        return (edgeMask/R).reshape((self.rows*self.columns,1))

    def getV(self, R):
        Cmat = self.getConductanceMatrix(R)
        Evec = self.getEdgeVec(R)
        return linsolve(Cmat, Evec)

    def Tdot(self, T, t):
        R = self.getResistance(T.reshape(self.rows,self.columns))
        V = self.getV(R)
        padded = np.pad(V, ((1, 1), (1, 1)), mode='constant')
        diffVsqr = (V - padded[2:, 1:-1])**2 + \
                   (V - padded[:-2, 1:-1])**2 + \
                   (V - padded[1:-1, 2:])**2 + \
                   (V - padded[1:-1, :-2])**2
        power = diffVsqr.flatten()/R.flatten()
        temperature_balance = self.tempConst.flatten()*(T**6 - self.bath_temperature**6)
        diffusion = self.diffusionConst*self.laplacian.dot(T.reshape((self.rows*self.columns,1)))
        return power - temperature_balance + diffusion.flatten()

    def getI(self, T):
        R = self.getResistance(T.reshape(self.rows,self.columns))
        V = self.getV(R)
        padded = np.pad(R, ((1, 1), (1, 1)), mode='constant')
        second_diagonal = 1 / (R + padded[1:-1, :-2])
        second_diagonal[:, 0] = 0
        second_diagonal = second_diagonal.flatten()
        second_diagonal = second_diagonal[1:]
        n_diagonal = 1 / (R[1:, :] + R[:-1, :])

        # Current coming from left
        diagonal = 1 / (R + padded[1:-1, :-2])
        matrix1 = -np.diagflat(diagonal) + \
                  np.diagflat(second_diagonal, k=-1)
        leftI = matrix1.dot(V.reshape((self.rows*self.columns,1))).reshape((self.rows, self.columns))
        leftI[:, 0] += self.VL/R[:, 0]
        # Current going to right
        diagonal = 1 / (R + padded[1:-1, 2:])
        matrix2 = np.diagflat(diagonal) - \
                 np.diagflat(second_diagonal, k=1)
        rightI = matrix2.dot(
            V.reshape((self.rows * self.columns, 1))).reshape(
            (self.rows, self.columns))
        rightI[:, -1] -= self.VR / R[:, -1]
        # Current coming from up
        diagonal = 1 / (R + padded[:-2, 1:-1])
        matrix3 = -np.diagflat(diagonal) + \
                 np.diagflat(n_diagonal, k=-self.columns)
        upI = matrix3.dot(
            V.reshape((self.rows * self.columns, 1))).reshape(
            (self.rows, self.columns))
        upI[0, :] += self.VU / R[0, :]
        # Current going down
        diagonal = 1 / (R + padded[2:, 1:-1])
        matrix4 = np.diagflat(diagonal) - \
                 np.diagflat(n_diagonal, k=self.columns)
        downI = matrix4.dot(
            V.reshape((self.rows * self.columns, 1))).reshape(
            (self.rows, self.columns))
        downI[-1, :] -= self.VD / R[-1, :]
        return leftI, rightI, upI, downI

    def setPlot(self, Imax, Imin, Tmax, Tmin):
        self.fig, self.ax = plt.subplots()
        self.horzIrows = np.arange(1, step=3, stop=self.rows*3)
        self.leftIcols = np.arange(0, step=3, stop=self.columns*3)
        self.rightIcols = np.arange(2, step=3, stop=self.columns * 3)
        self.upIrows = np.arange(0, step=3, stop=self.rows * 3)
        self.downIrows = np.arange(2, step=3, stop=self.rows * 3)
        self.vertIcols = np.arange(1, step=3, stop=self.columns * 3)
        self.Imask = np.ones((self.rows*3, self.columns*3))
        self.Imask[np.ix_(self.horzIrows, self.leftIcols)] = 0
        self.Imask[np.ix_(self.horzIrows, self.rightIcols)] = 0
        self.Imask[np.ix_(self.upIrows, self.vertIcols)] = 0
        self.Imask[np.ix_(self.downIrows, self.vertIcols)] = 0
        self.Tmask = np.ones((self.rows*3, self.columns*3))
        self.Tmask[np.ix_(self.horzIrows, self.vertIcols)] = 0

        Imasked = np.ma.masked_array(np.zeros((self.rows * 3, self.columns * 3)), self.Imask)
        self.Iimg = self.ax.imshow(Imasked, vmin=Imin, vmax=Imax)
        Tmasked = np.ma.masked_array(np.zeros((self.rows * 3, self.columns * 3)), self.Tmask)
        self.Timg = self.ax.imshow(Tmasked, vmin=Tmin, vmax=Tmax)

        self.text = self.ax.text(1, 1, "V = " + str(np.round(self.VL*100)/100))
        cb1 = plt.colorbar(self.Iimg, shrink=0.25)
        cb1.set_label('I')

        cb2 = plt.colorbar(self.Timg, shrink=0.25)
        cb2.set_label('T')

    def plotRes(self, result):
        V, I, T = result
        Ileft, Iright, Iup, Idown = I
        Iarray = np.zeros((3*self.rows, 3*self.columns))
        Iarray[np.ix_(self.horzIrows, self.leftIcols)] = Ileft
        Iarray[np.ix_(self.horzIrows, self.rightIcols)] = Iright
        Iarray[np.ix_(self.upIrows, self.vertIcols)] = Iup
        Iarray[np.ix_(self.downIrows, self.vertIcols)] = Idown
        I_masked = np.ma.masked_array(Iarray, self.Imask)
        self.Iimg.set_array(I_masked)

        T = T.reshape((self.rows, self.columns))
        Tarray = np.zeros((3*self.rows, 3*self.columns))
        Tarray[np.ix_(self.horzIrows, self.vertIcols)] = T
        T_masked = np.ma.masked_array(Tarray, self.Tmask)
        self.Timg.set_array(T_masked)

        self.text.set_text('V = ' + str(np.round(V*100)/100))

        return self.Iimg, self.Timg, self.text

    def plotICartoon(self, Vs, Is, Ts, path):
        '''
        updating the plot to current currents map
        :return: image for animation
        '''
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=24, bitrate=1800)
        self.setPlot(np.max(np.array(Is)),np.min(np.array(Is)),
                     np.max(np.array(Ts)),np.min(np.array(Ts)))
        frames = [(Vs[i], Is[i], Ts[i]) for i in range(len(Vs))]
        im_ani = animation.FuncAnimation(self.fig, self.plotRes,
                                         frames=frames, interval=100,
                                         repeat_delay=1000,
                                         blit=True)

        im_ani.save(path + '_cartoon.mp4', writer=writer)
        plt.close(self.fig)

    def plotIV(self, Vs, Is, path):
        Itot = []
        for I in Is:
            Ileft, Iright, Iup, Idown = I
            Itot.append(np.sum(Iright[:,-1]))
        mid = len(Vs)//2
        fig = plt.figure()
        plt.plot(Vs[:mid], Itot[:mid], '.b', label="increasing V")
        plt.plot(Vs[mid:], Itot[mid:], '.r', label="decreasing V")
        plt.legend()
        plt.xlabel("V")
        plt.ylabel("I")
        plt.savefig(path + "_IV.png")
        plt.close(fig)


    def getSteadyStateTemperature(self, T0):
        T = T0
        t = np.linspace(0, 1000, 10000)
        sol = odeint(array.Tdot, T, t)
        T = sol[-1]
        return T

    def getIV(self, vmax):
        steadyT = self.bath_temperature*np.ones((self.rows*self.columns,))
        Vvec = np.linspace(0,vmax,100)
        Vvec = np.hstack((Vvec,np.flipud(Vvec)))
        Is = []
        Ts = []
        for VL in Vvec:
            self.VL = VL
            steadyT = self.getSteadyStateTemperature(steadyT)
            I = self.getI(steadyT)
            plt.show()
            Is.append(I)
            Ts.append(steadyT)
            print(VL,end=',')
        return Vvec, Is, Ts

    def saveRes(self, Vs, Is, Ts):
        pass

    def saveParams(self):
        pass


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
    M = 30
    N = 30
    Tbath = 0.1
    navg = 10
    nstd = 9
    R0 = np.ones((M,N))
    R0[(0,-1),:] = 1000
    alpha = 0
    VL = 0
    VR = 0
    VU = 0
    VD = 0
    array = ThermalArray(M, N, navg, nstd, R0,
                         alpha, Tbath, VL, VR, VU, VD)
    dir_path = "dbg_thermal"
    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)
    file_name = "dbg_no_diffusion_30_30"
    path = os.path.join(dir_path,file_name)
    Vs, Is, Ts = array.getIV(3)
    array.plotICartoon(Vs, Is, Ts, path)
    array.plotIV(Vs, Is, path)

