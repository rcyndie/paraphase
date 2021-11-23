import numpy as np
import pytest
from paraphase.derivatives.j_compute import j_compute_slow
from paraphase.derivatives.jh_compute import jh_compute_slow
from paraphase.derivatives.jhr_compute import jhr_compute_slow
from paraphase.derivatives.jhr_compute import jhr_compute


def test_jhr():
    # defined coordinates
    nchan = 3
    freq = np.linspace(0.5, 1.5, nchan)
    ntime = 2
    time = np.linspace(0, 1, ntime)
    nsource = 3
    lcoord = 0.01 * np.random.randn(3)
    mcoord = 0.01 * np.random.randn(3)

    # freq profile
    def fprofile(v):
        return 1./v

    nant = 7
    ncorr = 2

    # parameters
    nparam = 3
    ntimeint = 1
    nchanint = 1
    alpha = 0.1 * np.random.randn(int(np.ceil(ntime/ntimeint)),
                                  int(np.ceil(nchan/nchanint)),
                                  nant, ncorr, nparam)

    def basis(l, m):
        return np.array([1, l, m])

    # # for non-parametric case
    # def basis(L, s):
    #     return L[s, :]

    # model vis
    model_arr = (np.random.randn(nsource, ntime, nchan, nant, nant, ncorr) +
                 1j*np.random.randn(nsource, ntime, nchan, nant, nant, ncorr))

    # vis
    vis_arr = np.zeros((ntime, nchan, nant, nant, ncorr), dtype=np.complex128)

    ##
    # residuals_arr = np.zeros((vis_arr.shape), dtype=vis_arr.dtype)

    for t in range(ntime):
        tt = t//ntimeint
        for f in range(nchan):
            ff = f//nchanint
            fprof = fprofile(freq[f])
            for p in range(nant):
                for q in range(nant):
                    for c in range(ncorr):
                        for s in range(nsource):
                            phasep = basis(lcoord[s], mcoord[s]).dot(alpha[tt, ff, p, c])
                            phaseq = basis(lcoord[s], mcoord[s]).dot(alpha[tt, ff, q, c])
                            vis_arr[t, f, p, q, c] += (np.exp(1j*fprof*phasep) *
                                    model_arr[s, t, f, p, q, c] * np.exp(-1j*fprof*phaseq))

    #for testing purposes only.
    residuals_arr = 0.5*vis_arr

    # evaluate jacobian with j_compute_slow and jh with jh_compute_slow.
    jacobian = j_compute_slow(vis_arr, model_arr, freq, alpha, fprofile, basis, lcoord, mcoord)
    jh = jh_compute_slow(vis_arr, model_arr, alpha, lcoord, mcoord, freq, fprofile, basis)

    # #Evaluate JHr explicitly with jhr_compute_slow().
    jhr = jhr_compute_slow(jh, residuals_arr, alpha)

    #Evaluate JHr.
    jhr_fast = jhr_compute(np.reshape(jh, (alpha.size, vis_arr.size)), np.reshape(residuals_arr, residuals_arr.size))
    jhr_fast = np.reshape(jhr_fast, alpha.shape)

    deltas = [1e-2, 1e-3]
    for tt in range(ntimeint):
        for ff in range(nchanint):
            for ant in range(nant):
                for cc in range(ncorr):
                    for par in range(nparam):                    
                        jhr2 = jhr[tt, ff, ant, cc, par]
                        jhr_fast2 = jhr_fast[tt, ff, ant, cc, par]
                        diffs = np.zeros(3)
                        for i, delta in enumerate(deltas):
                            diffs[i] = np.abs(jhr2 - jhr_fast2)
                        assert np.allclose(diffs[0]/diffs[1], deltas[0]/deltas[1],
                                        atol=1)


##Running tests.
test_jhr()