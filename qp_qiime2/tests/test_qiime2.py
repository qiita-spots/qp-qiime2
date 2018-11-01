# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os import remove
from shutil import rmtree
from tempfile import mkdtemp
from json import dumps
from os.path import exists, isdir, join, realpath, dirname

from qiita_client.testing import PluginTestCase

from qiime2 import __version__ as qiime2_version

from qp_qiime2 import plugin
from qp_qiime2 import call_qiime2


class qiime2Tests(PluginTestCase):
    def setUp(self):
        # this will allow us to see the full errors
        self.maxDiff = None

        plugin("https://localhost:21174", 'register', 'ignored')
        self._clean_up_files = []

        self.data = {
            'user': 'demo@microbio.me',
            'command': None,
            'status': 'running',
            'parameters': None}

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_rarefy(self):
        params = {
            'The feature table to be rarefied.': '5',
            'The total frequency that each sample should be rarefied to. '
            'Samples where the sum of frequencies is less than the sampling '
            'depth will be not be included in the resulting table.': '2',
            'qp-hide-method': u'rarefy',
            'qp-hide-paramThe total frequency that each sample should be '
            'rarefied to. Samples where the sum of frequencies is less than '
            'the sampling depth will be not be included in the resulting '
            'table.': 'sampling_depth',
            'qp-hide-paramThe feature table to be rarefied.': 'table',
            'qp-hide-plugin': 'feature-table'}
        self.data['command'] = dumps(
            ['qiime2', qiime2_version, 'Rarefy table'])
        self.data['parameters'] = dumps(params)

        jid = self.qclient.post(
            '/apitest/processing_job/', data=self.data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = call_qiime2(self.qclient, jid, params, out_dir)
        self.assertTrue(success)
        self.assertEqual(msg, '')
        self.assertEqual(ainfo[0].files, [(
            join(out_dir, 'rarefy', 'rarefied_table', 'feature-table.biom'),
            'biom')])
        self.assertEqual(ainfo[0].output_name, 'rarefied_table')

    def test_rarefy_error(self):
        params = {
            'The feature table to be rarefied.': '5',
            'The total frequency that each sample should be rarefied to. '
            'Samples where the sum of frequencies is less than the sampling '
            'depth will be not be included in the resulting table.': '200000',
            'qp-hide-method': u'rarefy',
            'qp-hide-paramThe total frequency that each sample should be '
            'rarefied to. Samples where the sum of frequencies is less than '
            'the sampling depth will be not be included in the resulting '
            'table.': 'sampling_depth',
            'qp-hide-paramThe feature table to be rarefied.': 'table',
            'qp-hide-plugin': 'feature-table'}
        self.data['command'] = dumps(
            ['qiime2', qiime2_version, 'Rarefy table'])
        self.data['parameters'] = dumps(params)

        jid = self.qclient.post(
            '/apitest/processing_job/', data=self.data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = call_qiime2(self.qclient, jid, params, out_dir)
        self.assertFalse(success)
        self.assertIsNone(ainfo)
        self.assertEqual(
            msg, 'Error: The rarefied table contains no samples or features. '
            'Verify your table is valid and that you provided a shallow '
            'enough sampling depth.')

    def test_beta(self):
        # now test phylogenetic
        params = {
            'Phylogenetic tree': join(
                dirname(realpath(__file__)), 'prune_97_gg_13_8.tre'),
            'The beta diversity metric to be computed.': 'unweighted_unifrac',
            'The feature table containing the samples over which beta '
            'diversity should be computed.': '5',
            '[Excluding weighted_unifrac] - The number of jobs to use for '
            'the computation. This works by breaking down the pairwise '
            'matrix into n_jobs even slices and computing them in parallel. '
            'If -1 all CPUs are used. If 1 is given, no parallel computing '
            'code is used at all, which is useful for debugging. For '
            'n_jobs below -1, (n_cpus + 1 + n_jobs) are used. Thus for '
            'n_jobs = -2, all CPUs but one are used. (Description from '
            'sklearn.metrics.pairwise_distances)': '1',
            'qp-hide-method': 'beta_phylogenetic',
            'qp-hide-paramPhylogenetic tree': 'phylogeny',
            'qp-hide-paramThe beta diversity metric to be computed.': 'metric',
            'qp-hide-paramThe feature table containing the samples over '
            'which beta diversity should be computed.': 'table',
            'qp-hide-plugin': 'diversity',
            'qp-hide_param[Excluding weighted_unifrac] - The number of jobs '
            'to use for the computation. This works by breaking down the '
            'pairwise matrix into n_jobs even slices and computing them '
            'in parallel. If -1 all CPUs are used. If 1 is given, no '
            'parallel computing code is used at all, which is useful for '
            'debugging. For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. '
            'Thus for n_jobs = -2, all CPUs but one are used. (Description '
            'from sklearn.metrics.pairwise_distances)': 'n_jobs'}
        self.data['command'] = dumps(
            ['qiime2', qiime2_version, 'Beta diversity (phylogenetic)'])
        self.data['parameters'] = dumps(params)

        jid = self.qclient.post(
            '/apitest/processing_job/', data=self.data)['job']
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = call_qiime2(self.qclient, jid, params, out_dir)
        self.assertEqual(msg, '')
        self.assertTrue(success)
        self.assertEqual(ainfo[0].files, [(
            join(out_dir, 'beta', 'distance_matrix', 'distance-matrix.tsv'),
            'plain_text')])
        self.assertEqual(ainfo[0].output_name, 'distance_matrix')

        # first let's test non phylogenetic
        params = {
            'A pseudocount to handle zeros for compositional metrics.  This '
            'is ignored for other metrics.': '5',
            'The beta diversity metric to be computed.': 'rogerstanimoto',
            'The feature table containing the samples over which beta '
            'diversity should be computed.': '8',
            'The number of jobs to use for the computation. This works '
            'by breaking down the pairwise matrix into n_jobs even slices '
            'and computing them in parallel. If -1 all CPUs are used. If '
            '1 is given, no parallel computing code is used at all, which '
            'is useful for debugging. For n_jobs below -1, (n_cpus + 1 + '
            'n_jobs) are used. Thus for n_jobs = -2, all CPUs but one are '
            'used. (Description from sklearn.metrics.pairwise_distances)': '1',
            'qp-hide-method': 'beta',
            'qp-hide-paramThe beta diversity metric to be computed.': 'metric',
            'qp-hide-paramThe feature table containing the samples over '
            'which beta diversity should be computed.': 'table',
            'qp-hide-plugin': 'diversity',
            'qp-hide_paramA pseudocount to handle zeros for compositional '
            'metrics.  This is ignored for other metrics.': 'pseudocount',
            'qp-hide_paramThe number of jobs to use for the computation. '
            'This works by breaking down the pairwise matrix into n_jobs even '
            'slices and computing them in parallel. If -1 all CPUs are used. '
            'If 1 is given, no parallel computing code is used at all, which '
            'is useful for debugging. For n_jobs below -1, (n_cpus + 1 + '
            'n_jobs) are used. Thus for n_jobs = -2, all CPUs but one are '
            'used. (Description from '
            'sklearn.metrics.pairwise_distances)': 'n_jobs'}
        self.data['command'] = dumps(
            ['qiime2', qiime2_version, 'Beta diversity'])
        self.data['parameters'] = dumps(params)

        jid = self.qclient.post(
            '/apitest/processing_job/', data=self.data)['job']
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = call_qiime2(self.qclient, jid, params, out_dir)
        self.assertEqual(msg, '')
        self.assertTrue(success)
        self.assertEqual(ainfo[0].files, [(
            join(out_dir, 'beta', 'distance_matrix', 'distance-matrix.tsv'),
            'plain_text')])
        self.assertEqual(ainfo[0].output_name, 'distance_matrix')


if __name__ == '__main__':
    main()
