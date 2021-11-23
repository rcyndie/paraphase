import numpy as np
from functools import partial
import pytest
from paraphase.derivatives.j_compute import j_compute_slow
from paraphase.derivatives.j_compute import j_compute

def test_jacobian(datadiag=True):
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

    # For nt_int = 1 and nf_int 1 we have
    # (nant - 1)/2 >= ndir
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

    # Test current implementation of Jacobian against j_compute_slow.
    jacobian2 = j_compute(vis_arr, model_arr, lcoord, mcoord, alpha, freq, fprofile, basis, datadiag)

    # test against finite differences
    deltas = [1e-2, 1e-3]
    for t in range(ntime):
        for f in range(nchan):
            fprof = fprofile(freq[f])
            for c in range(ncorr):
                for p in range(nant):
                    for q in range(nant):
                        # end of data axes
                        jac = jacobian[t, f, p, q, c]
                        jac2 = jacobian2[t, f, p, q, c]
                        model = model_arr[:, t, f, p, q, c]
                        # start of param axes
                        for tt in range(ntimeint):
                            for ff in range(nchanint):
                                for cc in range(ncorr):
                                    if cc==c and tt==t//ntimeint and ff==f//nchanint and p!=q:
                                        # test against finite difference approx here
                                        vis_func = partial(vis_model, p=p, q=q,
                                                           alphas=alpha[tt, ff, :, cc], model=model,
                                                           lcoord=lcoord, mcoord=mcoord, fprof=fprof,
                                                           basis=basis)
                                        for ant in range(nant):
                                            if (ant==p or ant==q):
                                                for par in range(nparam):
                                                    j = jac[tt, ff, ant, cc, par]
                                                    j2 = jac2[tt, ff, ant, cc, par]
                                                    derivs = np.zeros(3, dtype=model.dtype)
                                                    diffs = np.zeros(3)
                                                    diffs2 = np.zeros(3)
                                                    for i, delta in enumerate(deltas):
                                                        vishigh = vis_func(ant, par, 1, delta)
                                                        vislow = vis_func(ant, par, -1, delta)
                                                        deriv = (vishigh - vislow)/(2*delta)
                                                        diffs[i] = np.abs(j - 2*deriv)
                                                        diffs2[i] = np.abs(j2 - j)
                                                    # assert np.allclose(diffs[0]/diffs[1], deltas[0]/deltas[1],
                                                                    #    atol=1)
                                                    assert np.allclose(diffs2[0]/diffs2[1], deltas[0]/deltas[1],
                                                                       atol=1)


def vis_model(ant, par, sign, delta, p, q, alphas, model, lcoord, mcoord, fprof, basis):
    # print(ant, par, sign)
    nsource = model.size
    vis = 0j
    alpha = alphas[ant]
    alpha[par] += sign * delta
    for s in range(nsource):
        if ant==p:
            L = basis(lcoord[s], mcoord[s])
            gp = np.exp(1j * fprof * L.dot(alpha))
            gq = np.exp(1j * fprof * L.dot(alphas[q]))
            vis += gp * model[s] * np.conj(gq)
        elif ant==q:
            L = basis(lcoord[s], mcoord[s])
            gp = np.exp(1j * fprof * L.dot(alphas[p]))
            gq = np.exp(1j * fprof * L.dot(alpha))
            vis += gp * model[s] * np.conj(gq)
    return vis


##Running tests.
test_jacobian()
