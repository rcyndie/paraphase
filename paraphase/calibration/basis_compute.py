import numpy as np

def _make_basis_vec(n_param, l_s, m_s):
    """
    Generating the basis polynomial to compute the phase equation. Right now, 
    it is of the form [1, l_s, m_s, l_s**2, m_s**2, ...] with l_s and m_s being 
    scalars. The function returns a vector of length n_params.
    """
    
    N = (n_param+1) // 2
    lvec = (np.tile(l_s, N-1) ** np.arange(1, N))
    mvec = (np.tile(m_s, N-1) ** np.arange(1, N))

    main_vec = np.ones((N-1, 2))
    main_vec[:, 0] = lvec
    main_vec[:, 1] = mvec
    main_vec = (main_vec).flatten()
    
    return np.insert(main_vec, 0, [1])

# def _get_basis(n_param, sources):
#     """
#     Get basis matrix of shape (n_params, n_dir). Both l and m, which represent the
#     direction cosine coordinates of the sources are each vectors of length n_dir.
    
#     """

#     #Get the dimension and the direction cosines.
#     l = sources[:, 1]
#     m = sources[:, 2]

#     #Get "n_dir" dimension.
#     n_dir = sources.shape[0] #len(l)
#     basis = np.zeros((n_param, n_dir))

#     for s in range(n_dir):
#         basis[:, s] = _make_basis_vec(n_param, l[s], m[s])

#     return basis

def _get_basis(sources):
    """
    Get basis matrix >>> covariance matrix >>> Cholesky decomposition.
    
    """

    #Get covariance matrix.
    K = _squared_exp(sources, sources)

    #
    jitter = 1e-6

    #Compute Cholesky decomposition of K_inv.
    L = np.linalg.cholesky(K + jitter*np.eye(K.shape[0]))
    #Remember L and K are of shapes n_sources \times n_sources, or
    #if I am incorrect, they should at least be of the same shapes.

    return L

def _squared_exp(x, xp, sigmaf=1000, l=0.1):
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