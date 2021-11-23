import numpy as np
from paraphase.format_1D2D import get_xxyy

def residuals_compute(data_arr, model_arr, lcoord, mcoord, alpha, chan_freq, freqf, basis, datadiag=False):
    """
    The function returns difference between model_arr and data_arr
    at this instance.

    """

    #Initialise residuals.
    residuals = data_arr.copy()
    
    n_timint, n_freint, n_ant, n_ccor, _ = alpha.shape
    n_tim = data_arr.shape[0]
    n_fre = data_arr.shape[1]
    n_dir = model_arr.shape[0]

    t_int = n_tim//n_timint
    f_int = n_fre//n_freint
    
    for t in range(n_tim):
        tt = t//t_int
        for f in range(n_fre):
            ff = f//f_int
            fprof = freqf(chan_freq[f])
            for p in range(n_ant):
                for q in range(p): #note only doing this for q < p
                    for d in range(n_dir):
                        for k in range(n_ccor):
                            #Subtract model for each direction.
                            L = basis(lcoord[d], mcoord[d])
                            gp = np.exp(1j * fprof * L.dot(alpha[tt, ff, p, k]))
                            gq = np.exp(1j * fprof * L.dot(alpha[tt, ff, q, k]))
                            if datadiag:
                                residuals[t, f, p, q, k] -= gp * model_arr[d, t, f, p, q, k] * np.conj(gq)
                            else:
                                residuals[t, f, p, q, k, k] -= gp * model_arr[d, t, f, p, q, k, k] * np.conj(gq)

                    residuals[t, f, q, p] = np.conj(residuals[t, f, p, q])

    return np.reshape(residuals, residuals.size)