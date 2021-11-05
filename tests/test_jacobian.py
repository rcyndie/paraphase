import numpy as np
from functools import partial
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
    def basis(l, m):
        return np.array([1, l, m])

    # model vis
    ncorr = 1
    model_arr = (np.random.randn(ntime, nchan, nant, nant, ncorr, nsource) +
                 1j*np.random.randn(ntime, nchan, nant, nant, ncorr, nsource))

    # vis
    vis_arr = np.zeros((ntime, nchan, nant, nant, ncorr), dtype=np.complex128)
    for t in range(ntime):
        for f in range(nchan):
            fprof = fprofile(freq[f])
            for p in range(nant):
                for q in range(nant):
                    for c in range(ncorr):
                        for s in range(nsource):
                            phasep = basis(lcoord[s], mcoord[s]).dot(alpha[t, f, p, c])
                            phaseq = basis(lcoord[s], mcoord[s]).dot(alpha[t, f, q, c])
                            vis_arr[t, f, p, q, c] += (np.exp(1j*fprof*phasep) *
                                    model[t, f, p, q, c, s] * np.exp(-1j*fprof*phaseq))

    # evaluate jacobian
    jacobian = j_compute_slow(data_arr, model_arr, freq, alpha, fprofile, basis, lcoord, mcoord)

    # test against finite differences
    delta_param = 1e-3
    for t in range(ntime):
        for f in range(nchan):
            fprof = fprofile(freq[f])
            for c in range(ncorr):
                for p in range(nant):
                    for q in range(nant):
                        # end of data axes
                        jac = jacobian[t, f, p, q, c]
                        model = model_arr[t, f, p, q, c]
                        # start of param axes
                        for tt in range(ntimeint):
                            for ff in range(nchanint):
                                for cc in range(ncorr):
                                    # test against finite difference approx here
                                    vis_func = partial(vis_model, p=p, q=q,
                                                        alphas=alpha[tt, ff, :, cc], model=model,
                                                        lcoord=lcoord, mcoord=mcoord, fprof=fprof,
                                                        basis=basis, delta=delta)
                                    for ant in range(nant):
                                        for par in range(npar):
                                            vislow = vis_func(ant, par, 1)
                                            vishigh = vis_func(ant, par, -1)
                                            deriv = (vishigh - vislow)/(2*delta_param)
                                            j = jac[tt, ff, ant, cc, par]
                                            print(np.abs(j - deriv))
                                            assert np.allclose(j, deriv, atol=delta_alpha**2)




def vis_model(ant, par, sign, p, q, alphas, model, lcoord, mcoord, fprof, basis, delta):
    nsource = model.size
    vis = 0j
    alpha = alphas[ant]
    for s in range(nsource):
        alpha = alphas[ant]
        alpha[par] += sign * delta
        if ant == p:
            L = basis(lcoord[s], mcoord[s])
            Lalpha = L.dot(alpha)
            prefact = basis(lcoord[s], mcoord[s])[par]
            phase = 1j * fprof * Lalpha
            gp = np.exp(phase)
            gq = np.exp(1j*fprof*L.dot(alphas[q]))
            vis += prefact * phase * gp * model[s] * np.conj(gq)
        elif ant == q:
            L = basis(lcoord[s], mcoord[s])
            Lalpha = L.dot(alpha)
            prefact = basis(lcoord[s], mcoord[s])[par]
            phase = 1j * fprof * Lalpha
            gq = np.exp(phase)
            gp = np.exp(1j*fprof*L.dot(alphas[p]))
            vis += -prefact * phase * gp * model[s] * np.conj(gq)
    return vis
