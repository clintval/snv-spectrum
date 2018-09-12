import setuptools

from pathlib import Path
from setuptools import find_packages

PACKAGE = 'nucleic'
VERSION = '0.5.0'

setuptools.setup(
    name='nucleic',
    version=VERSION,
    author='clintval',
    author_email='valentine.clint@gmail.com',
    description='A Python 3.6 library for plotting mutational spectra.',
    url=f'https://github.com/clintval/{PACKAGE}',
    download_url=f'https://github.com/{PACKAGE}/archive/v{VERSION}.tar.gz',
    long_description=Path('README.md').read_text(),
    long_description_content_type='text/markdown',
    license='MIT',
    zip_safe=False,
    packages=[PACKAGE],
    install_requires=['biopython', 'matplotlib', 'mpl_helpers', 'numpy', 'palettable', 'pyfaidx'],
    keywords='signature mutation transition transversion spectra bioinformatics',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    project_urls={
        'Documentation': f'https://{PACKAGE}.readthedocs.io',
        'Issue-Tracker': f'https://github.com/clintval/{PACKAGE}/issues',
    },
)
