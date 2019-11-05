# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main

from qiime2.sdk import PluginManager

from qp_qiime2.util import get_qiime2_type_name_and_predicate


class UtilTests(TestCase):
    def test_get_qiime2_type_name_and_predicate(self):
        pm = PluginManager()
        q2plugin = pm.plugins['feature-classifier']
        parameters = q2plugin.methods['classify_sklearn'].signature.parameters

        # testing regular parameters
        exp = ('Int', None)
        obs = get_qiime2_type_name_and_predicate(parameters['n_jobs'])
        self.assertEqual(exp, obs)

        # testing union paramters
        exp = ('Float', {'type': 'predicate', 'name': 'Range',
                         'range': [0, 1], 'inclusive': [True, True]})
        obs = get_qiime2_type_name_and_predicate(parameters['confidence'])
        self.assertEqual(exp, obs)


if __name__ == '__main__':
    main()
