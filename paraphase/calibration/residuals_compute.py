import numpy as np

def residuals_compute(data_arr, model_arr, gains, diag=False):
    """
    The function returns difference between model_arr and data_arr
    at this instance.

    """

    #Initialise residuals.
    residuals = data_arr.copy()

    n_tim, n_fre, n_ant, _, n_ccor, _ = data_arr.shape
    n_dir = gains.shape[0]
    n_timint = gains.shape[1]

    for t in range(n_tim):
        tt = t//n_timint
        for f in range(n_fre):
            for p in range(n_ant):
                for q in range(p): #note only doing this for q < p
                    for d in range(n_dir):
                        #Subtract model for each direction.
                        residuals[t, f, p, q] -= gains[d, tt, f, p] * model_arr[d, t, f, p, q] * np.conj(gains[d, tt, f, q].T)

    if diag is True:
        return get_xx_yy_residual(residuals)
    else:
        return np.reshape(residuals, (n_tim*n_fre*n_ant*n_ant*n_ccor*n_ccor))

def get_xx_yy_residual(residuals):
    """
    Extract only the XX and YY components from the residuals.

    """

    #
    n_tim, n_fre, n_ant, _, n_ccor, _ = residuals.shape

    #Initialise the new residuals.
    new_residual = np.zeros((n_tim, n_fre, n_ant, n_ant, n_ccor), dtype=residual.dtype)

    for k in range(n_ccor):
        new_residual[..., k] = residual[..., k, k]

    return new_residual.reshape(n_tim*n_fre*n_ant*n_ant*n_ccor)