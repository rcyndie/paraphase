import numpy as np
import pytest
from paraphase.derivatives.j_compute import j_compute_slow
from paraphase.derivatives.jh_compute import jh_compute
from paraphase.derivatives.jh_compute import jh_compute_slow

def test_jh():
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

    # evaluate jacobian with j_compute_slow.
    jacobian = j_compute_slow(vis_arr, model_arr, freq, alpha, fprofile, basis, lcoord, mcoord)

    # Test current implementation of complex conjugate transpose of Jacobian against jh_compute_slow.
    jh = jh_compute_slow(vis_arr, model_arr, alpha, lcoord, mcoord, freq, fprofile, basis)
    jh_fast = jh_compute(np.reshape(jacobian, (vis_arr.size, alpha.size)))
    jh_fast = np.reshape(jh_fast, alpha.shape+vis_arr.shape)

    deltas = [1e-2, 1e-3]
    for t in range(ntime):
        for f in range(nchan):
            fprof = fprofile(freq[f])
            for c in range(ncorr):
                for p in range(nant):
                    for q in range(nant):
                        # end of data axes
                        jac = jacobian[t, f, p, q, c]
                        jh2 = jh[:, :, :, :, :, t, f, p, q, c]
                        jh_fast2 = jh_fast[:, :, :, :, :, t, f, p, q, c]
                        model = model_arr[:,t, f, p, q, c]
                        # start of param axes
                        for tt in range(ntimeint):
                            for ff in range(nchanint):
                                for cc in range(ncorr):
                                    if cc==c and tt==t//ntimeint and ff==f//nchanint and p!=q:
                                        for ant in range(nant):
                                            if (ant==p or ant==q):
                                                for par in range(nparam):
                                                    jach = np.conjugate(jac[tt, ff, ant, cc, par])
                                                    jh3 = jh2[tt, ff, ant, cc, par]
                                                    jh_fast3 = jh_fast2[tt, ff, ant, cc, par]
                                                    diffs = np.zeros(3)
                                                    diffs2 = np.zeros(3)
                                                    for i, delta in enumerate(deltas):
                                                        diffs[i] = np.abs(jach - jh3)
                                                        diffs2[i] = np.abs(jh3 - jh_fast3)
                                                    assert np.allclose(diffs[0]/diffs[1], deltas[0]/deltas[1],
                                                                    atol=1)
                                                    assert np.allclose(diffs2[0]/diffs2[1], deltas[0]/deltas[1],
                                                                    atol=1)


##Running tests.
test_jh()
