import numpy as np
from pyrap.tables import table
from pyrap.tables import maketabdesc
from pyrap.tables import makearrcoldesc

def ms_1D_to_2D(msname, column, tchunk=None, fchunk=1, n_dir=1, DD=False):
    """
	The function reshapes the visibilites from an input ms into a shape which
    can be passed to a CubiCal solver
    (n_mod, n_time, n_fre, n_ant, n_ant, n_cor, n_cor).

    """

    def _convert_to_2D(data, N):
        """
		convert the 1D timeslots data from the MS to 2D NxN data

		"""

        data2x2 = np.reshape(data,(data.shape[0], 2, 2))
        arr = np.zeros((N, N, 2, 2), dtype=data.dtype)
        indices = np.triu_indices(N, 1)
        rind, cind = np.tril_indices(N, -1)
        arr[indices] = data2x2
        arr[rind, cind] = arr[cind, rind].conj()
        return arr

    tab = table(msname)
    n_ant = max(tab.getcol("ANTENNA2", 0, -1, 1)) + 1
    nbl = n_ant*(n_ant-1)/2

    NTimes = len(np.unique(tab.getcol("TIME")))
    if tchunk is None:
        tchunk = NTimes

    ntchunks = NTimes//tchunk

    if DD is False:
        vis = {}
        vis["src"+str(0)] = tab.getcol(column)
        ntot, nfre, ncc = vis["src"+str(0)].shape
    else:
        vis = {}
        for i in range(n_dir):
            vis["src"+str(i)] = tab.getcol(column+str(i))
            ntot, nfre, ncc = vis["src"+str(0)].shape

    tab.close()

    temparray = np.zeros((NTimes, n_dir, nfre, n_ant, n_ant, 2, 2), dtype=vis["src0"].dtype)
    for d in range(n_dir):
        for i in xrange(NTimes):
            for f in range(nfre):
                temparray[i,d,f,:,:,:,:] = _convert_to_2D(vis["src"+str(d)][i*nbl:(i+1)*nbl,f,:], n_ant)

    if DD:
        outarray = np.zeros((ntchunks, n_dir, tchunk, nfre, n_ant, n_ant, 2, 2), dtype=temparray.dtype)
        for d in range(n_dir):
            outarray[:,d,...] = np.reshape(temparray[:,d,...],(ntchunks, tchunk, nfre, n_ant, n_ant, 2, 2))
    else:
        outarray = np.reshape(temparray,(ntchunks, tchunk, nfre, n_ant, n_ant, 2, 2))

    return outarray


def ms_2D_to_1D(msname, column, in_array, tchunk=120, fchunk=1, chan=False, timerow=False, valuetype=None):
    """
    Reshapes back the data MS style and store the output in the corresponding MS column.
    
    """

    #Add the output column if it doesn't exist.
    valuetype = valuetype or 'complex'
    # im.argo.addcol(msname=msname, colname=column, valuetype=valuetype)
    add_col(msname, colname=column, valuetype=valuetype, initialiseto=0., clone="DATA")

    def _convert_to_1D(data,N):
        """
        convert the 2D timeslots data to 1D data

        """

        data4 = np.reshape(data, (data.shape[0], data.shape[1], 4))
        arr = data4[np.triu_indices(N, 1)]
        return arr

    tab = table(msname, readonly=False)
    print(tab.colnames())
    n_ant = max(tab.getcol("ANTENNA2")) + 1
    nbl = n_ant*(n_ant-1)/2
    vis = tab.getcol(column) 

    NTimes = len(np.unique(tab.getcol("TIME")))

    ntchunks, ntim, nfre, n_ant, n_ant, n_cor, n_cor = in_array.shape
    ntt = ntchunks*ntim

    newarray = np.reshape(in_array,(ntt,nfre,n_ant,n_ant,2,2))

    outarray = np.zeros((ntt*nbl, nfre, 4), dtype=vis.dtype)
    for i in xrange(ntt):
        for f in xrange(nfre):
            outarray[i*nbl:(i+1)*nbl,f,:] = _convert_to_1D(newarray[i,f,:,:,:,:], n_ant)

    if chan:
        vis[:,0:nfre,:] = outarray[:,0:nfre,:]
    elif timerow:
        vis[0:ntt*nbl,:,:] = outarray
    else:
        vis = outarray

    tab.putcol(column, vis)
    tab.close()


def add_col(msname, colname=None, shape=None, valuetype=None, initialiseto=None, coldmi=None, rowchunk=None, clone="DATA"):
    """
    The function adds a new column colname to the MS.

    msname - Name of MS
    colname - Name of column to be created or rewritten
    shape - Shape of column
    valuetype - Type of the values in the colname
    initialiseto - Initialise colname to this value
    coldmi - dminfo of existing column (getdminfo())
    rowchunk - Step size or size of chunk
    clone - Alternative column to be cloned from
    
    """

    ms = table(msname, readonly=False)

    if colname in ms.colnames():
        print("The column %s already exists!"%colname)
    
    #getcell(columnname, rownr) << get data from a column cell.
    if shape:
        data_desc = maketabdesc(makearrcoldesc(colname, value=initialiseto, shape=shape, ndim=0, valuetype=valuetype))
    elif clone:
        element = ms.getcell(clone, 0)
        shape = element.shape
        data_desc = maketabdesc(makearrcoldesc(colname, value=0., shape=shape, ndim=0, valuetype=valuetype))

    colinfo = [data_desc, coldmi] if coldmi else [data_desc]
    ms.addcols(*colinfo)
    print("The column %s has been successfully added."%colname) 
    
    if initialiseto is None:
        ms.close()
    else:
        spwids = set(ms.getcol("DATA_DESC_ID"))
        for spw in spwids:
            ms_spw = ms.query("DATA_DESC_ID={0:d}".format(spw))
            nrows = ms_spw.nrows()
            rowchunk = nrows/10
            dshape = [0] + [a for a in shape]
            for row0 in range(0, nrows, rowchunk):
                nr = min(rowchunk, nrows-row0)
                dshape[0] = nr
                ms_spw.putcol(colname, np.ones(dshape, dtype=type(initialiseto))*initialiseto, row0, nr)
            ms_spw.close()

    ms.close()


def get_xxyy(arr1):
    """
    Extract only the autocorrelations from the array arr1.

    """
    
    #
    n_ccor = arr1.shape[-1]

    #Initialise the new array.
    arr2 = np.zeros(arr1.shape[:-1], dtype=arr1.dtype)

    for k in range(n_ccor):
        arr2[..., k] = arr1[..., k, k]

    return arr2