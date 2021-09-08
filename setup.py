from setuptools import setup, find_packages
#No clue about find_packages as of now

#import paraphase

with open("README.md", "r") as f:
	long_description = f.read()

requirements = [
				'python-casacore',
				'argparse',
				]



setup(
	name='paraphase',
	version="0.0.1",
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
	package_dir={"": "paraphase"},
	#packages=["paraphase"],
	entry_points='''
					[console_scripts]
					paraphase=paraphase.paraphase:main
	'''
	,
)
