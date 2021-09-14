import numpy as np

def gains_compute():
    """
    Returns the gains of shape (n_dir, n_ant) using alpha (array containing gain
    parameters per antenna).
    """ 

    #Initialise gains.
    #self.gains = np.zeros((self.n_dir, self.n_timint, self.n_fre, self.n_ant, self.n_cor, self.n_cor), dtype=self.dtype)

    # for s in range(self.n_dir):
    #     for t in range(self.n_timint):
    #         for f in range(self.n_fre):
    #             ff = f//self.f_int
    #             for p in range(self.n_ant):
    #                 for k in range(self.n_cor):
    #                     #To correct the dimension issue ((n_param,) instead of (n_param, 1)).
    #                     alpha_vec = (self.alpha[t, ff, p, :, k]).reshape(self.n_param)
    #                     #phase_equation = np.dot(alpha_vec, basis[:, s])
    #                     self.gains[s, t, f, p, k, k] = np.exp(1.0j * self.chunk_fs[f] * np.dot(alpha_vec, self.basis[:, s]))

    for s in range(self.n_dir):
        for t in range(self.n_timint):
            for f in range(self.n_fre):
                ff = f//self.f_int
                for p in range(self.n_ant):
                    for k in range(self.n_cor):
                        self.gains[s, t, f, p, k, k] = np.exp(1.0j * self.chunk_fs[f] *  self.basis[s].dot(self.alpha[t, ff, p, :, k]))


