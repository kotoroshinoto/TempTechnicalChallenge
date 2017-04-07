import sys
import os
from setuptools import setup, find_packages
from pip.req import parse_requirements

if sys.version_info.major < 3:
	print("I'm only for python 3, please upgrade")
	sys.exit(1)

install_reqs = parse_requirements(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'requirements.txt'), session=False)
reqs = [str(ir.req) for ir in install_reqs]

setup(
	name='annotate_vcf_challenge',
	version='0.0.1',
	packages=find_packages(),
	include_package_data=True,
	install_requires=reqs,
	description="annotate variants from a VCF file using INFO values and ExAC database",
	long_description="""\
	programming challenge, annotating VCF file
	""",
	author="Michael Gooch",
	author_email="goochmi@gmail.com",
	url="https://github.com/kotoroshinoto/TempTechnicalChallenge",
	classifiers=[
		"Programming Language :: Python :: 3 :: Only",
	],
	entry_points={
		'console_scripts': [
		'annotate_vcf_challenge = annotate_vcf_challenge.annotate_vcf:annotate'
		]
	}
)
