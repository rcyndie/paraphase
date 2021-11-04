import pytest
from paraphase.derivatives.j_compute import j_compute_slow

def test_jacobian():
    # defined coordinates
    nchan = 3
    freq = np.linspace(0.5, 1.5, nchan)
    ntime = 2
    time = np.linspace(0, 1, ntime)
    nsource = 3
    lcoord = 0.1 * np.random.randn(3)
    mcoord = 0.1 * np.random.randn(3)

    # freq profile
    def fprofile(v):
        return 1./v

    # For nt_int = 1 and nf_int 1 we have
    # (nant - 1)/2 >= ndir
    nant = 7

    # parameters
    nparam = 3
    alpha = np.random.randn(ntime, nfreq, nant, ncorr, nparam)
    def basis(l, m, alp):
        return 1*alp[0] + l*alp[1] + m*alp[2]

    # model vis
    ncorr = 1
    model_arr = (np.random.randn(ntime, nchan, nant, nant, ncorr, nsource) +
                 1j*np.random.randn(ntime, nchan, nant, nant, ncorr, nsource))

    # vis
    vis_arr = np.zeros((ntime, nchan, nant, nant, ncorr), dtype=np.complex128)
    for t in range(ntime):
        for f in range(nchan):
            for p in range(nant):
                for q in range(nant):
                    for c in range(ncorr):
                        fprof = fprofile(freq[f])
                        for s in range(nsource):
                            phasep = basis(lcoord[s], mcoord[s], alpha[t, f, p, c])
                            phaseq = basis(lcoord[s], mcoord[s], alpha[t, f, q, c])
                            vis_arr[t, f, p, q, c] += (np.exp(1j*fprof*phasep) *
                                                       model_arr[t, f, p, q, c, s] *
                                                       np.exp(-1j*fprof*phaseq))

    # evaluate jacobian
    jac = j_compute_slow(data_arr, model_arr, freq, alpha, fprofile, basis, l, m)

    # test against finite differences
