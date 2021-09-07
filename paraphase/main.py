from __future__ import print_function
import argparse
import yaml
from yaml.loader import SafeLoader
import sys
# from paraphase.paraphase.calibration.calibrate import ParametrisedPhase


def create_parser(argv=None):
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
	#about gains
	p1 = p.add_argument_group("g")
	p1.add_argument("--g-save-to", dest="save-to", type=str, help="Save gains to this address")
	p1.add_argument("--g-time-int", dest="time-int", type=int, help="Size of solution time interval")
	p1.add_argument("--g-freq-int", dest="freq-int", type=int, help="Size of solution frequency interval")
	#about data
	p2 = p.add_argument_group("data")
	p2.add_argument("--data-ms", dest="ms", help="Name of measurement set", type=str, action="append")
	p2.add_argument("--data-column", dest="column", type=str, help="Name of MS column to read for data")
	#about sky
	p3 = p.add_argument_group("sky")
	p3.add_argument("--sky-model", dest="sky-model", type=str, help="Tigger lsm file", action="append")
	#about output
	p4 = p.add_argument_group("out")
	p4.add_argument("--out-writeto", dest="out-writeto", type=str, help="Write to output MS column")
	
	args = p.parse_args(remaining_argv)

	return args


def main(debugging=False):
	"""
	Main paraphase driver function.

	"""

	args = create_parser()
	print(args)
	


