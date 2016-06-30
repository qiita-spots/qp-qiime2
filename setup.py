# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
from setuptools import setup, find_packages

setup(
    name="qp_qiime2",
    version="0.0.0-dev",
    packages=find_packages(),
    install_requires=['q2-feature-table', 'scikit-bio',
                      'qiime >= 2.0.0', 'click', 'qiita_client', 'h5py'],
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="A QIIME 2 plugin for Qiita",
    license="BSD",
    scripts=glob.glob('scripts/*'),
    url="http://www.qiime.org",
)
