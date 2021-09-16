import numpy as np
import configparser
import argparse
# import pyyaml
# from yaml.loader import SafeLoader
import sys
import paraphase
from optparse import OptionParser
from pyrap.tables import table
from paraphase.calibration.calibrate import calibratewith

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
	p.add_option("--gtype", dest="gtype", type=str, help="Specify basis solver for parametrised phase gains")
	p.add_option("--deltachi", dest="deltachi", type=float, help="Specify threshold for solution stagnancy")
	# p.add_option("--msname", dest="msname", help="Name of measurement set", action="append")
	p.add_option("--column", dest="column", type=str, help="Name of MS column to read for data")
	p.add_option("--model", dest="model", type=str, help="Tigger lsm file", action="append")
	p.add_option("--writeto", dest="writeto", type=str, help="Write to output MS column")
	
	return p

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

	data = options.column
	msrcs = options.model

	if len(args) != 1:
		ri('Please specify a Measurement Set to calibrate.')
		sys.exit(-1)
	else:
		#Remove any trailing characters, for example "/".
		msname = args[0].rstrip('/')


	#MS info.
	# spwtab = table(msname+"/FIELD")
	# n_chan = spwtab.getcol("NUM_CHAN")
	# spwtab.close()

	freqtab = table(msname+"/SPECTRAL_WINDOW")
	chan_freq = freqtab.getcol("CHAN_FREQ").squeeze()
	freqtab.close()
	n_fre = chan_freq.size

	anttab = table(msname+"/ANTENNA")
	n_ant = len(anttab)
	antnames = anttab.getcol("NAME")
	anttab.close()

	tt = table(msname)
	uniants = np.unique(tt.getcol("ANTENNA1"))
	tt.close()

	###
	#maybe write the gains shape and alpha shape.
	n_timint = options.timint
	n_freint = options.freint
	n_dir = len(msrcs)
	print(n_dir, "length of the model list")

	##How to get the following?
	n_cor = 2
	n_param = 3
	# n_dir = 3

	solvertype = options.gtype
	alpha_shape = [n_timint, n_freint, n_ant, n_param, n_cor]
	gains_shape = [n_dir, n_timint, n_fre, n_ant, n_cor, n_cor]
	calibratewith(data, msrcs, options.deltachi, alpha_shape, gains_shape, options.deltachi, solvertype)


	return options, args


if __name__ == "__main__":
	options, args = main()
	# print(args[0])
	# print(options[0])
