from setuptools import setup, find_packages
#No clue about find_packages as of now

import paraphase

with open("README.md", "r") as f:
	long_description = f.read()

requirements = [
				'python-casacore',
				'argparse',
<<<<<<< HEAD
				'yaml',
=======
				'pyyaml',
>>>>>>> 2c47f081839836d3a77a382b6d29a89454b2537d
				]



setup(
	name='paraphase',
	version=paraphase.__version__,
	author="Cyndie Russeeawon",
	author_email="crusseeaw@gmail.com",
	description="Parametrised phase direction-dependent calibration",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/rcyndie/paraphase",
	project_urls={
		"Bug Tracker": "https://github.com/rcyndie/paraphase/issues",
	},
	install_requires=requirements,
	classifiers=[
			"Programming Language :: Python :: 3",
			"License :: OSI Approved :: MIT License",
			"Operating System :: OS Independent",
	],
<<<<<<< HEAD
	#package_dir={"": "paraphase"},
	#packages=["paraphase"],
	entry_points='''
					[console_scripts]
					paraphase=paraphase.paraphase.main:main
=======
	# package_dir={"": "paraphase"},
	packages=find_packages(),
	entry_points='''
					[console_scripts]
					paraphase=paraphase.main:main
>>>>>>> 2c47f081839836d3a77a382b6d29a89454b2537d
	'''
	,
)
