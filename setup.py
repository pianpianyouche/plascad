#!/usr/bin/env python3.6
from setuptools import setup
setup(name='Plascad',
      version='1.16',
      packages=['plas_cad'],
      description="Plasmid classfication, ARGs annotation, Visualization",
      author="You CHE",
      author_email="cheyou@hku.hk",
      url='https://github.com/pianpianyouche/plascad',
      keywords=['Plasmids', 'antibiotic resistance', 'mobilization'],
      license='MIT',
      install_requires=['biopython'],
      include_package_data=True,
      package_dir={'plascad': 'plas_cad'},
      package_data={'plas_cad': ['scripts/*','scripts/*/*','database/*/*','database/*/*/*','example/*']},
      entry_points={'console_scripts': ['Plascad = plas_cad.__main__:main']},
      classifiers=['Programming Language :: Python :: 3']
      )
