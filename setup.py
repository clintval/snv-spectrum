import io

from setuptools import find_packages, setup

with io.open('README.md', encoding='utf-8') as f:
    long_description = f.read()

__version__ = '0.3.0'

DOWNLOAD_URL = 'https://github.com/clintval/snv-spectrum/archive/v{}.tar.gz'

KEYWORDS = [
    'bioinformatics',
    'mutation',
    'signature',
    'spectra',
    'transition',
    'transversion']

setup(
    name='snv_spectrum',
    packages=find_packages(),
    version=__version__,
    description='A Python 3.6 library for plotting mutational spectra.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='clintval',
    author_email='valentine.clint@gmail.com',
    url='https://github.com/clintval/snv-spectrum',
    download_url=DOWNLOAD_URL.format(__version__),
    install_requires=[
        'biopython',
        'matplotlib',
        'mpl_helpers',
        'numpy',
        'palettable',
        'pyfaidx'
    ],
    extras_require={
        'ci': ['nose', 'codecov'],
        'cluster-ready': ['fastcluster', 'scipy', 'polo'],
    },
    license='MIT',
    zip_safe=True,
    keywords=(
        'signature mutation transition transversion spectra bioinformatics'
    ),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.6',
    ]
)
