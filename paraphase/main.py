import logging
import numpy as np
import configparser
import argparse
# import pyyaml
# from yaml.loader import SafeLoader
import sys
import paraphase
import Tigger
from optparse import OptionParser
from pyrap.tables import table
from paraphase.calibration.calibrate import calibratewith
from paraphase.format_1D2D import ms_1D_to_2D
from paraphase.format_1D2D import ms_2D_to_1D
from paraphase.format_1D2D import get_xxyy

#Create a custom logger.
# logger = logging.getLogger(__name__)
# logger.warning("This is a warning")

def ri(message):
	"""
	'\033[91m' << display text in red.

	"""

	print('\033[91m'+message+'\033[0m')


def create_parser0(argv=None):
	"""
	Create an ArgumentParser object.

	"""

	if argv is None:
		argv = sys.argv

	pconf = argparse.ArgumentParser(description=__doc__,
								formatter_class=argparse.RawDescriptionHelpFormatter,
								add_help=False)

	pconf.add_argument("-c", "--conf_file", help="Specify configuration file", metavar="FILE", default="default.yml")
	args, remaining_argv = pconf.parse_known_args()

	default = {}

	if args.conf_file:
		# c = ConfigParser.SafeConfigParser()
		config = configparser.ConfigParser()
		config.read([args.conf_file])
		default.update(dict(config.items("data")))
		default.update(dict(config.items("g")))
		default.update(dict(config.items("sky")))
		default.update(dict(config.items("out")))


	#Parsing the rest of the arguments.
	p = argparse.ArgumentParser(parents=[pconf])

	p.set_defaults(**default)
	#about data
	p1 = p.add_argument_group("data")
	p1.add_argument("--data-ms", dest="ms", help="Name of measurement set", type=str, action="append")
	p1.add_argument("--data-column", dest="column", type=str, help="Name of MS column to read for data")
	#about gains
	p2 = p.add_argument_group("g")
	p2.add_argument("--g-save-to", dest="save-to", type=str, help="Save gains to this address")
	p2.add_argument("--g-type", type=str, help="Specify basis for parametrised phase gains")
	p2.add_argument("--g-time-int", dest="time-int", type=int, help="Size of solution time interval")
	p2.add_argument("--g-freq-int", dest="freq-int", type=int, help="Size of solution frequency interval")
	#about sky
	p3 = p.add_argument_group("sky")
	p3.add_argument("--sky-model", dest="sky-model", type=str, help="Tigger lsm file", action="append")
	#about output
	p4 = p.add_argument_group("out")
	p4.add_argument("--out-writeto", dest="out-writeto", type=str, help="Write to output MS column")

	args = p.parse_args(remaining_argv)

	return args


def create_parser(argv=None):
	"""
	Create a parser object with OptionParser().

	"""

	p = OptionParser(usage='%prog [options] msname')
	# p.add_option("-c", "--config", help="Specify configuration file", metavar="FILE", default="default.yml", action="append")
	p.add_option("-l", "--list", dest="dolist", action="store_true", help="List MS properties and exit", default=False)
	p.add_option("--save-to", dest="save-to", type=str, help="Save gains to this address")
	p.add_option("--timint", dest="timint", type=int, help="Size of solution time interval")
	p.add_option("--freint", dest="freint", type=int, help="Size of solution frequency interval")
	p.add_option("--freqf", dest="freqf", type=str, help="Specify function for frequency dependence", default="linear")
	p.add_option("--gtype", dest="gtype", type=str, help="Specify basis solver for parametrised phase gains", default="ppoly")
	p.add_option("--npar", dest="npar", type=int, help="Specify (odd) number of parameters for gtype-ppoly", default=3)
	p.add_option("--kernel", dest="kernel", type=str, help="Specify kernel for covariance function for gtype-pcov")
	p.add_option("--sigmaf", dest="sigmaf", type=float, help="Standard deviation which controls vertical scaling for gtype-pcov", default=1000.0)
	p.add_option("--lscale", dest="lscale", type=float, help="Specify input length-scale for gtype-pcov", default=1.0)
	p.add_option("--jitter", dest="jitter", type=float, help="Specify jitter for Cholesky decomposition for gtype-pcov", default=1e-6)
	p.add_option("--deltachi", dest="deltachi", type=float, help="Specify threshold for solution stagnancy")
	# p.add_option("--msname", dest="msname", help="Name of measurement set", action="append")
	p.add_option("--column", dest="column", type=str, help="Name of MS column to read for data", default="DATA")
	p.add_option("--datadiag", dest="datadiag", action="store_true", help="Specify for diagonal data", default=False)
	p.add_option("--skymodel", dest="skymodel", type=str, help="Tigger lsm file")
	p.add_option("--writeto", dest="writeto", type=str, help="Write to output MS column")
	
	return p

def extract_modelsrcs(phase_centre, model):
	"""
	How to manipulate Tigger sky model (lsm) (model) to obtain ra, dec?
	>>> l and m.

	"""

	mdsrcs = Tigger.load(model)
	n_dir = len(mdsrcs)
	ra0, dec0 = phase_centre[0], phase_centre[1]
	arr1 = np.zeros((n_dir, 2), dtype=np.float64)
	arr2 = np.zeros((n_dir, 3), dtype=np.float64)

	for d in range(n_dir):
		arr2[d, 0] = mdsrcs[d].flux.I
		arr1[d, 0] = mdsrcs[d].pos.ra
		arr1[d, 1] = mdsrcs[d].pos.dec

	ra_delta = arr1[:, 0] - ra0
	arr2[:, 1] = np.cos(arr1[:, 1]) * np.sin(ra_delta)
	arr2[:, 2] = np.sin(arr1[:, 1]) * np.cos(dec0) - np.cos(arr1[:, 1]) * np.sin(dec0) * np.cos(ra_delta)

	return arr2


def debug():
	"""

	"""

	main(debugging=True)


def main(debugging=False):
	"""
	Main paraphase driver function.

	"""

	#Create parser object.
	(options, args) = create_parser().parse_args()

	skymodel = options.skymodel

	if len(args) != 1:
		ri('Please specify a Measurement Set to calibrate.')
		sys.exit(-1)
	else:
		#Remove any trailing characters, for example "/".
		msname = args[0].rstrip('/')

	#MS info.
	fieldtab = table(msname+"/FIELD")
	phase_centre = fieldtab.getcol("PHASE_DIR")[0, 0]
	fieldtab.close()

	freqtab = table(msname+"/SPECTRAL_WINDOW")
	chan_freq = freqtab.getcol("CHAN_FREQ").squeeze()
	freqtab.close()
	n_fre = chan_freq.size
	if n_fre > 1:
		chan_freq = chan_freq/min(chan_freq)
		if options.freqf == "linear":
			chan_freq = 1./chan_freq
		elif options.freqf == "quadratic":
			chan_freq = 1/(chan_freq**2)
	elif n_fre == 1:
		chan_freq = np.ones(1, dtype=chan_freq.dtype)

	anttab = table(msname+"/ANTENNA")
	n_ant = len(anttab)
	antnames = anttab.getcol("NAME")
	anttab.close()

	tt = table(msname)
	uniants = np.unique(tt.getcol("ANTENNA1"))
	#Calibrate *** column.
	data = tt.getcol(options.column)
	tt.close()
	n_cor = data.shape[2]
	#Consider a 2 \times 2 shape for n_cor.
	#Let n_ccor * n_ccor = n_cor.
	n_ccor = n_cor//2

	#Specify solution interval sizes.
	n_timint = options.timint
	n_freint = options.freint

	#About modelled sources.
	ra0, dec0 = phase_centre[0], phase_centre[1]
	arr_srcs = extract_modelsrcs(phase_centre, skymodel)
	n_dir = arr_srcs.shape[0]

	#Parameters for basis.
	bparams = {"gtype": options.gtype}
	if options.gtype == "ppoly":
		n_par = options.npar
		bparams["n_par"] = options.npar
	elif options.gtype == "pcov":
		n_par = n_dir
		bparams2 = {"n_par": n_dir, "sigmaf": options.sigmaf, "lscale": options.lscale, "jitter": options.jitter, "kernel": options.kernel}
		bparams.update(bparams2)
	
	#Parameters for gains.
	#Consider diagonal gains.
	alpha_shape = [n_timint, n_freint, n_ant, n_par, n_ccor]
	gains_shape = [n_dir, n_timint, n_fre, n_ant, n_ccor]
	gparams = {"chan_freq": chan_freq, "alpha_shape": alpha_shape, "gains_shape": gains_shape}
	gparams.update(bparams)

	#Parameters for solutions.
	sparams = {"deltachi": options.deltachi}

	#
	data_arr = ms_1D_to_2D(msname, column="DATA", tchunk=1, fchunk=1, n_dir=1, DD=False)
	data_arr = data_arr[:, 0]
	model_arr = ms_1D_to_2D(msname, "DD_SRC_", tchunk=None, fchunk=1, n_dir=n_dir, DD=True)
	model_arr = model_arr[0]

	if options.datadiag:
		data_arr = get_xxyy(data_arr)
		model_arr = get_xxyy(model_arr)
	
	# ms_2D_to_1D(msname, column="DATA3", in_array=data, tchunk=1, fchunk=1, chan=False, timerow=False, valuetype=None)

	calibratewith(data_arr, model_arr, arr_srcs, bparams, gparams, sparams, options.datadiag)

	return options, args


if __name__ == "__main__":
	options, args = main()