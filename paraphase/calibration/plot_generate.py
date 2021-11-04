import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def chi2_plot(arr1, sparams):
    """
    Plotting arr1.
    
    """

    #Fontsize.
    a = 20
    matplotlib.rcParams.update({'font.size': a})
    
    niter = len(arr1)
    fig = plt.figure(figsize=(15, 8))
    plt.plot(arr1)
    plt.xlabel("Iteration number")
    plt.xticks(np.arange(1, niter, 1))
    plt.ylabel(r"$\chi^2$")
    plt.tight_layout()
    fig.savefig(sparams["outputdir"]+"/Chi2.pdf")
    plt.close("all")

def jhj_plot(jhj, sparams):
    """
    Plotting JHJ.

    """
    
    fig = plt.figure('jhj real')
    plt.title("Re(JHJ)")
    plt.imshow(jhj.real)
    plt.colorbar()
    fig.savefig(sparams["outputdir"]+"/jhj_re.pdf")
    plt.close("all")

def jac_plot(jac, sparams):
    """
    Plotting Jacobian.
    
    """

    fig = plt.figure('Jacobian')
    plt.title("Re(J)")
    plt.imshow(jac.real, aspect="auto")
    plt.colorbar()
    fig.savefig(sparams["outputdir"]+"/jac_re.pdf")
    plt.close("all")

    fig = plt.figure('Jacobian')
    plt.title("Im(J)")
    plt.imshow(jac.imag, aspect="auto")
    plt.colorbar()
    fig.savefig(sparams["outputdir"]+"/jac_im.pdf")
    plt.close("all")
