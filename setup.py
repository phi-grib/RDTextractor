from setuptools import setup
import glob, sys

data_files_list = glob.glob('data/*pkl')
data_files_list.extend(glob.glob('data/*pkl.gz'))

data_files = [
    ('data/', data_files_list)
]

setup(name='rdt_extractor',
      version='0.1',
      license='GNU GPLv3',
      description='Script to filter and extract data from the eTox Vitic database.',
      url='https://github.com/phi-grib/RDTextractor',
      download_url='https://github.com/phi-grib/RDTextractor.git',
      author='Elisabet Gregori and Ignacio Pasamontes',
      author_email='elisabet.gregori@upf.edu',
      packages=['rdt_extractor'],
      scripts=['rdt_extractor/extract'],
      data_files=data_files
    )
