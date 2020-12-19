#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
# Standard imports
#
import glob, os
#
# setuptools' sdist command ignores MANIFEST.in
#
from setuptools import setup, find_packages
#
#
# Begin setup
#
setup_keywords = dict()
#
# THESE SETTINGS NEED TO BE CHANGED FOR EVERY PRODUCT.
#
setup_keywords['name'] = 'astropath'
setup_keywords['description'] = 'Probabilistic Association of Transients to Hosts (PATH)'
setup_keywords['author'] = 'FRB Associations Team'
setup_keywords['author_email'] = 'xavier@ucolick.org'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/FRBs/astropath'

#
# END OF SETTINGS THAT NEED TO BE CHANGED.
#
setup_keywords['version'] = '0.1.dev0' #get_version(setup_keywords['name'])
#
# Use README.rst as long_description.
#
setup_keywords['long_description'] = ''
if os.path.exists('README.md'):
    with open('README.md') as readme:
        setup_keywords['long_description'] = readme.read()
#
# Set other keywords for the setup function.  These are automated, & should
# be left alone unless you are an expert.
#
# Treat everything in bin/ except *.rst as a script to be installed.
#
if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in glob.glob(os.path.join('bin', '*'))
        if not os.path.basename(fname).endswith('.rst')]
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['requires'] = ['Python (>3.7.0)']
setup_keywords['zip_safe'] = False
setup_keywords['packages'] = find_packages()
#setup_keywords['package_dir'] = {'':'py'}
setup_keywords['setup_requires']=['pytest-runner']
setup_keywords['tests_require']=['pytest']

#
# Add internal data directories.
#

data_files = []

# walk through the data directory, adding all files
data_generator = os.walk('astropath/data')
for path, directories, files in data_generator:
    for f in files:
        data_path = '/'.join(path.split('/')[1:])
        data_files.append(data_path + '/' + f)
setup_keywords['package_data'] = {'astropath': data_files,
                                  '': ['*.rst', '*.txt', '*.yaml']}
setup_keywords['include_package_data'] = True

#
# Run setup command.
#
setup(**setup_keywords)

