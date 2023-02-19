# Author: Yudai Okubo <yudaiokubo@gmail.com>
# Copyright (c) 2022-2023 Yudai Okubo
# License: MIT

from setuptools import setup

DESCRIPTION = "Inter-atomic quantity reporter for OpenMM"
NAME = 'curp-reporter'
AUTHOR = 'Yudai Okubo'
AUTHOR_EMAIL = 'yudaiokubo@gmail.com'
URL = 'https://github.com/passive-radio/openmm-curp-reporter'
LICENSE = 'MIT'
DOWNLOAD_URL = 'https://github.com/passive-radio/openmm-curp-reporter'
VERSION = '0.0.1'
PYTHON_REQUIRES = '>=3.6'
KEYWORDS = 'openmm curp md'

INSTALL_REQUIRES = [
    'numpy==1.24.1',
    'OpenMM==8.0.0'
]

PACKAGES = [
    'reporter', 'test', 'input'
]

CLASSIFIERS = [
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
]

with open('README.md', 'r', encoding='utf-8') as fp:
    readme = fp.read()

LONG_DESCRIPTION = readme
LONG_DESCRIPTION_CONTENT_TYPE = 'text/markdown'

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type=LONG_DESCRIPTION_CONTENT_TYPE,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    maintainer=AUTHOR,
    maintainer_email=AUTHOR_EMAIL,
    url=URL,
    download_url=URL,
    packages=PACKAGES,
    classifiers=CLASSIFIERS,
    license=LICENSE,
    keywords=KEYWORDS,
    include_package_data=True,
    install_requires=INSTALL_REQUIRES
)