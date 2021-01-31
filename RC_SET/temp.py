import numpy as np
import os

if __name__ == "__main__":
    directory = "2d_long_array_bgu_with_perp_point_contact_it"
    for run in [1, 2, 3, 4]:
        Imap = np.load(os.path.join(directory,"array_5_15_r_std_9_run_%d_it_Imap.npy" % run))
        I = []
        vertI = []
        for frame in Imap:
            i = (np.sum(frame[1::2, 0]) + np.sum(frame[1::2, -1])) / 2
            i_vert = (np.sum(frame[0, :]) + np.sum(frame[-1, :])) / 2
            I.append(i)
            vertI.append(np.abs(i_vert))
        np.save(os.path.join(directory,'array_5_15_r_std_9_run_%d_it_I.npy' % run), I)
        np.save(os.path.join(directory,'array_5_15_r_std_9_run_%d_it_vertI.npy' % run), vertI)