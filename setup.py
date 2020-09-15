#!/usr/bin/env python3
from setuptools import setup
setup(name='plas-cad',
      version='1.0',
      description="Plasmid classfication, ARGs annotation, Visualization",
      author="You CHE",
      author_email="cheyou@hku.hk",
      url='https://github.com/pianpianyouche/plas-cad',
      keywords=['Plasmids', 'antibiotic resistance', 'mobilization'],
      license='MIT',
      install_requires=['biopython'],
      python_requires='>=3.6',
      include_package_data=True,
      package_dir={'plas_cad': 'plas_cad'},
      package_data={'plas_cad': ['scripts/*','scripts/*/*','database/*/*','database/*/*/*','example/*']},
      entry_points={'console_scripts': ['plas-cad = plas_cad.__main__:main']},
      classifiers=['Programming Language :: Python :: 3']
      )

