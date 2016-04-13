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
    install_requires=['feature_table', 'scikit-bio >= 0.4.2, < 0.5.0',
                      'qiime >= 2.0.0', 'requests', 'click'],
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="A QIIME 2 plugin for Qiita",
    license="BSD",
    scripts=glob.glob('scripts/*'),
    url="http://www.qiime.org",
)
