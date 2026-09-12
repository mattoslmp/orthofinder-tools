# based on https://realpython.com/pypi-publish-python-package
# How to upload:
#  - change package version in `setup.py` and `__init__.py`
#  - `python setup.py sdist`
#  - `twine upload dist/orthofinder-tools-?.tar.gz`
import pathlib
from setuptools import setup

HERE = pathlib.Path(__file__).parent
README = (HERE / 'README.md').read_text(encoding='utf-8')

setup(
    name='orthofinder-tools',
    version='0.1.0',
    description='Annotate OrthoFinder orthogroups and create publication-ready comparative-genomics figures',
    long_description=README,
    long_description_content_type='text/markdown',
    url='https://github.com/mattoslmp/orthofinder-tools/',
    author='Thomas Roder; Databiomics extensions',
    license='MIT',
    packages=['orthofinder_tools'],
    python_requires='>=3.9',
    install_requires=[
        'numpy',
        'pandas',
        'biopython',
        'fire',
        'matplotlib',
        'seaborn'
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
    entry_points={
        'console_scripts': [
            'annotate_orthogroups=orthofinder_tools.annotate_orthogroups:main',
            'orthofinder_plots=orthofinder_tools.orthofinder_plots:main',
        ],
    },
)
