import numpy as np

def get_gains_poly(basis, alpha, gparams):
    """
    Returns the gains of shape (n_dir, n_ant) using alpha (array containing gain
    parameters per antenna).
    #Diagonal gains.
    
    """ 

    n_dir, n_timint, n_fre, n_ant, n_ccor = gparams["gains_shape"]
    n_par = gparams["n_par"]
    n_freint = gparams["alpha_shape"][1]
    gains = np.empty(gparams["gains_shape"], dtype=complex)

    for d in range(n_dir):
        for t in range(n_timint):
            for f in range(n_fre):
                ff = f//n_freint
                for p in range(n_ant):
                    for k in range(n_ccor):
                        gains[d, t, f, p, k] = np.exp(1.0j * gparams["chan_freq"][f] * np.dot(alpha[t, ff, p, :, k], basis[:, d]))

    return gains

def get_gains_cov(basis, alpha, gparams):
    """
    Returns the gains of shape (n_dir, n_ant) using alpha (array containing gain
    parameters per antenna).
    #Diagonal gains.
    
    """ 

    n_dir, n_timint, n_fre, n_ant, n_ccor = gparams["gains_shape"]
    n_freint = gparams["alpha_shape"][1]
    gains = np.empty(gparams["gains_shape"], dtype=complex)

    for d in range(n_dir):
        for t in range(n_timint):
            for f in range(n_fre):
                ff = f//n_freint
                for p in range(n_ant):
                    for k in range(n_ccor):
                        gains[d, t, f, p, k] = np.exp(1.0j * gparams["chan_freq"][f] * basis[d].dot(alpha[t, ff, p, :, k]))

    return gains

def gains_compute(basis, alpha, gparams):
    """
    The function computes a basis given the specifications.
    params is n_par when using a gtype-ppoly.
    params is n_dir, sigmaf, l and more when using a gtype-pcov.

    """

    if gparams["gtype"] == "ppoly":
        return get_gains_poly(basis, alpha, gparams)
    elif gparams["gtype"] == "pcov":
        return get_gains_cov(basis, alpha, gparams)