from setuptools import setup


try:
    import pypandoc
    long_description = pypandoc.convert_file('README.md', 'rst')
    long_description = long_description.replace('\r', '')
except (ImportError, OSError):
    import io
    with io.open('README.md', encoding='utf-8') as f:
        long_description = f.read()


setup(
    name='snv_spectrum',
    packages=['snv_spectrum'],
    version='0.2.0',
    description='A Python 3.6 library for plotting mutational spectra.',
    long_description=long_description,
    author='clintval',
    author_email='valentine.clint@gmail.com',
    url='https://github.com/clintval/snv-spectrum',
    download_url='https://github.com/clintval/snv-spectrum/archive/v0.2.0.tar.gz',
    py_modules=['snv_spectrum'],
    install_requires=[
        'biopython',
        'matplotlib',
        'numpy',
    ],
    extras_require={
        'ci': ['nose', 'codecov'],
        'fancytest': ['nose', 'nose-progressive', 'coverage'],
    },
    license='MIT',
    zip_safe=True,
    keywords='signature mutation transition transversion spectra bioinformatics',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.6',
    ]
)
