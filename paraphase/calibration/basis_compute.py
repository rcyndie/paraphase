import numpy as np

def make_basis_vec(n_par, l_d, m_d):
    """
    Generating the basis polynomial to compute the phase equation. Right now, 
    it is of the form [1, l_d, m_d, l_d**2, m_d**2, ...] with l_d and m_d being 
    scalars. The function returns a vector of length n_par.

    """
    
    N = (n_par+1) // 2
    lvec = (np.tile(l_d, N-1) ** np.arange(1, N))
    mvec = (np.tile(m_d, N-1) ** np.arange(1, N))

    main_vec = np.ones((N-1, 2))
    main_vec[:, 0] = lvec
    main_vec[:, 1] = mvec
    main_vec = (main_vec).flatten()
    
    return np.insert(main_vec, 0, [1])

def get_basis_poly(sources, n_par):
    """
    Get basis matrix of shape (n_params, n_dir). Both l and m, which represent the
    direction cosine coordinates of the sources are each vectors of length n_dir.
    
    """

    #Get the dimension and the direction cosines.
    l = sources[:, 1]
    m = sources[:, 2]

    n_dir = sources.shape[0]
    basis = np.zeros((n_par, n_dir))

    for d in range(n_dir):
        basis[:, d] = make_basis_vec(n_par, l[d], m[d])

    return basis

def get_basis_cov(sources, bparams):
    """
    Get basis matrix >>> covariance matrix >>> Cholesky decomposition.
    recall bparams = {"n_par": n_dir, "sigmaf": options.sigmaf, "lscale": options.lscale}
    and more.

    """

    #Get covariance matrix.
    #(sources, sources)? (x, xp)? maybe needs tweaking?
    K = squared_exp(sources, sources, bparams["sigmaf"], bparams["lscale"])

    #Compute Cholesky decomposition of K_inv.
    L = np.linalg.cholesky(K + bparams["jitter"] * np.eye(K.shape[0]))
    #Remember L and K are of shapes n_sources \times n_sources, or
    #if I am incorrect, they should at least be of the same shapes.

    return L

def squared_exp(x, xp, sigmaf, l):
    """
    The function returns a covariance matrix using the squared 
    exponential (SE) kernel.
    
    k(x, xp) = sigma^2 * exp(-(x-xp)^2/(2*l^2))
    
    """
    
    # if x.ndim > 1 or xp.ndim > 1:
    #    raise ValueError("Inputs must be 1D")
    
    #Get shape of x and xp.
    N = x.shape[0]
    M = xp.shape[0]

    #Create covariance matrix.
    C = np.zeros((N, M))

    for i in range(N):
        for j in range(M):
            C[i, j] = sigmaf**2*np.exp(-(1/2*l**2)*((x[:, 1][i] - xp[:, 1][j])**2 + (x[:, 2][i] - xp[:, 2][j])**2))
            
    return C


def basis_compute(msrcs, bparams):
    """
    The function computes a basis given the specifications.
    params is n_par when using a gtype-ppoly.
    params is n_dir, sigmaf, l and more when using a gtype-pcov.

    """

    if bparams["gtype"] == "ppoly":
        return get_basis_poly(msrcs, bparams["n_par"])
    elif bparams["gtype"] == "pcov":
        return get_basis_cov(msrcs, bparams)