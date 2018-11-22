import setuptools
import sys

from pathlib import Path
from setuptools import find_packages

PACKAGE: str = 'nucleic'

if sys.version_info < (3, 6):
    sys.exit(f'Python < 3.6 will not be supported.')


def this_version() -> str:
    """Read the variable `__version__` from the module itself."""
    contents = Path('nucleic/_version.py').read_text()
    *_, version = contents.strip().split()
    return version


setuptools.setup(
    name=PACKAGE,
    version=this_version(),
    author='clintval',
    author_email='valentine.clint@gmail.com',
    description='A Python 3.6 library for plotting mutational spectra.',
    url=f'https://github.com/clintval/{PACKAGE}',
    download_url=f'https://github.com/clintval/{PACKAGE}/archive/v{this_version()}.tar.gz',
    long_description=Path('README.md').read_text(),
    long_description_content_type='text/markdown',
    license='MIT',
    packages=setuptools.find_packages(where='./'),
    install_requires=[
        'fastcluster',
        'nimfa',
        'matplotlib',
        'numpy',
        'ordered-set',
        'palettable',
        'polo',
        'pyfaidx',
        # TODO: In skbio>0.5.0, NumPy is a build dependency
        'scikit-bio>=0.4.0,<=0.5.0',
        'scipy>=1.0.0',
        'toyplot',
    ],
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
        'Issue Tracker': f'https://github.com/clintval/{PACKAGE}/issues',
    },
    zip_safe=False,
)
