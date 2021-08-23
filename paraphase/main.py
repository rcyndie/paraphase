from __future__ import print_function
import configparser
import yaml
from yaml.loader import SafeLoader

def create_parser(argv=None):
	"""
	Create a ArgumentParser object.

	"""

	# p0 = argparse.ArgumentParser(
	# 		description=__doc__,
	# 		formatter_class=argparse.RawDescriptionHelpFormatter,
	# 		add_help=False)
	p0 = argparse.ArgumentParser()
	
	p = p0.add_argument_group("Input")
	p.add_argument("--ms", dest="ms", help="Name of measurement set", type=str, action="append")
	p.add_argument("--sky_model", dest="sky_model", type=str, help="Tigger lsm file", action="append")
	p.add_argument("--cfg", dest="cfg", type=argparse.FileType(mode='r'))

	return p0


def main(debugging=False):
	"""
	Main paraphase driver function.

	"""

	args = create_parser().parse_args()
	


