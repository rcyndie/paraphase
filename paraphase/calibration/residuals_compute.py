import numpy as np
from paraphase.format_1D2D import get_xxyy

def residuals_compute(data_arr, model_arr, gains, datadiag=False):
    """
    The function returns difference between model_arr and data_arr
    at this instance.

    """

    #Initialise residuals.
    residuals = data_arr.copy()

    n_dir, n_timint, n_fre, n_ant, n_ccor = gains.shape
    n_timint = gains.shape[1]
    n_tim = data_arr.shape[0]

    for t in range(n_tim):
        tt = t//n_timint
        for f in range(n_fre):
            for p in range(n_ant):
                for q in range(p): #note only doing this for q < p
                    for d in range(n_dir):
                        for k in range(n_ccor):
                            #Subtract model for each direction.
                            if datadiag:
                                residuals[t, f, p, q, k] -= gains[d, tt, f, p, k] * model_arr[d, t, f, p, q, k] * np.conj(gains[d, tt, f, q, k].T)
                            else:
                                residuals[t, f, p, q, k, k] -= gains[d, tt, f, p, k] * model_arr[d, t, f, p, q, k, k] * np.conj(gains[d, tt, f, q, k].T)

                    residuals[t, f, q, p] = np.conj(residuals[t, f, p, q])

    return np.reshape(residuals, residuals.size)