import argparse
import configparser
import os.path
import sys
import yaml
from yaml.loader import SafeLoader
import json


class LoadFromFile(argparse.Action):
    def __call__ (self, parser, namespace, values, option_string = None):
        with values as f:
            # parse arguments in the file and store them in the target namespace
            data = parser.parse_args(f.read().split(), namespace=None)
        for k, v in vars(data).items():
        	#set arguments in the target namespace if they have not been set yet,
        	if getattr(namespace, k, None) is not None:
        		setattr(namespace, k, v)


def is_valid_file(parser, arg):
	"""
	Check if the file exists in specified directory.

	"""

	if not os.path.exists(arg):
		parser.error("The file %s does not exist."%arg)
	else:
		return open(arg, 'r')


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


# #Is this redundant??
def parse_args(parser):
	"""
	Load the configuration file fed on the command line,
	and parse.

	"""	

	args = parser.parse_args()
	if args.cfg:
		try:
			data = yaml.load(args.cfg, Loader=SafeLoader)
		except yaml.YAMLError as exc:
			print(exc)

		delattr(args, 'cfg')
		arg_dict = args.__dict__
		for key, value in data.items():
			if isinstance(value, list):
				for v in value:
					arg_dict[key].append(v)
			else:
				arg_dict[key] = value

	return args


def create_ms(args):
	"""
	Create an empty ms with given arguments. To be filled later!

	"""

	os.system("simms -T kat-7 -t ascii -n {} -st 0.5 -dt 30 -ra 11h49m36s -dec -30d16m41s --nchan 5 -f0 1.4GHz -df 20MHz --pol 'XX XY YX YY' /home/russeeawon/Cyndie_jake/simms/simms/observatories/kat-7.itrf.txt".format(...))


if __name__ == "__main__":
	args = parse_args(create_parser())
	#Then, feed 'args' in some function.

	print(args.ms)
	#print(args.ms+args.sky_model)
	#And, for example,
	#create_ms(args)

