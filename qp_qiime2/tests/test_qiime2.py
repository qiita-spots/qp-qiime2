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

    def test_not_analysis_artifact(self):
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
        self.assertFalse(success)
        self.assertEqual(msg, 'Artifact "5" is not an analysis artifact.')

    def test_rarefy(self):
        params = {
            'The feature table to be rarefied.': '8',
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
            'The feature table to be rarefied.': '8',
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
            msg, 'Error running: The rarefied table contains no samples or '
            'features. Verify your table is valid and that you provided a '
            'shallow enough sampling depth.')

    def test_beta(self):
        params = {
            'A pseudocount to handle zeros for compositional metrics.  This '
            'is ignored for other metrics.': '1',
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
            'qp-hide-paramA pseudocount to handle zeros for compositional '
            'metrics.  This is ignored for other metrics.': 'pseudocount',
            'qp-hide-paramThe number of jobs to use for the computation. '
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

    def test_beta_phylogenetic(self):
        params = {
            'Phylogenetic tree': join(
                dirname(realpath(__file__)), 'prune_97_gg_13_8.tre'),
            'The beta diversity metric to be computed.': 'unweighted_unifrac',
            'The feature table containing the samples over which beta '
            'diversity should be computed.': '8',
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
            'qp-hide-param[Excluding weighted_unifrac] - The number of jobs '
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
            join(out_dir, 'beta_phylogenetic', 'distance_matrix',
                 'distance-matrix.tsv'), 'plain_text')])
        self.assertEqual(ainfo[0].output_name, 'distance_matrix')

    def test_alpha(self):
        params = {
            'The alpha diversity metric to be computed.': 'simpson',
            'The feature table containing the samples for which alpha '
            'diversity should be computed.': '8',
            'qp-hide-method': 'alpha',
            'qp-hide-paramThe alpha diversity metric to be '
            'computed.': 'metric',
            'qp-hide-paramThe feature table containing the samples for '
            'which alpha diversity should be computed.': 'table',
            'qp-hide-plugin': 'diversity'}
        self.data['command'] = dumps(
            ['qiime2', qiime2_version, 'Alpha diversity'])
        self.data['parameters'] = dumps(params)

        jid = self.qclient.post(
            '/apitest/processing_job/', data=self.data)['job']
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = call_qiime2(self.qclient, jid, params, out_dir)
        self.assertEqual(msg, '')
        self.assertTrue(success)
        self.assertEqual(ainfo[0].files, [(
            join(out_dir, 'alpha', 'alpha_diversity', 'alpha-diversity.tsv'),
            'plain_text')])
        self.assertEqual(ainfo[0].output_name, 'alpha_diversity')

    def test_alpha_phylogenetic(self):
        params = {
            'Phylogenetic tree': join(
                dirname(realpath(__file__)), 'prune_97_gg_13_8.tre'),
            'The alpha diversity metric to be computed.': 'faith_pd',
            'The feature table containing the samples for which alpha '
            'diversity should be computed.': '8',
            'qp-hide-method': 'alpha_phylogenetic',
            'qp-hide-paramPhylogenetic tree': 'phylogeny',
            'qp-hide-paramThe alpha diversity metric to be '
            'computed.': 'metric',
            'qp-hide-paramThe feature table containing the samples for '
            'which alpha diversity should be computed.': 'table',
            'qp-hide-plugin': 'diversity'}
        self.data['command'] = dumps(
            ['qiime2', qiime2_version, 'Alpha diversity (phylogenetic)'])
        self.data['parameters'] = dumps(params)

        jid = self.qclient.post(
            '/apitest/processing_job/', data=self.data)['job']
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = call_qiime2(self.qclient, jid, params, out_dir)
        self.assertEqual(msg, '')
        self.assertTrue(success)
        self.assertEqual(ainfo[0].files, [(
            join(out_dir, 'alpha_phylogenetic', 'alpha_diversity',
                 'alpha-diversity.tsv'), 'plain_text')])
        self.assertEqual(ainfo[0].output_name, 'alpha_diversity')

    def test_alpha_correlation(self):
        # as we don't have an alpha vector available, we will calculate
        # one using a non phylogenetic metric
        params = {
            'The alpha diversity metric to be computed.': 'simpson',
            'The feature table containing the samples for which alpha '
            'diversity should be computed.': '8',
            'qp-hide-method': 'alpha',
            'qp-hide-paramThe alpha diversity metric to be '
            'computed.': 'metric',
            'qp-hide-paramThe feature table containing the samples for '
            'which alpha diversity should be computed.': 'table',
            'qp-hide-plugin': 'diversity'}
        self.data['command'] = dumps(
            ['qiime2', qiime2_version, 'Alpha diversity'])
        self.data['parameters'] = dumps(params)

        jid = self.qclient.post(
            '/apitest/processing_job/', data=self.data)['job']
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = call_qiime2(self.qclient, jid, params, out_dir)
        data = {'filepaths': dumps(ainfo[0].files), 'type': "alpha_vector",
                'name': "Non phylogenetic alpha vector", 'analysis': 1,
                'data_type': '16S'}
        reply = self.qclient.post('/apitest/artifact/', data=data)
        aid = reply['artifact']

        params = {
            'The correlation test to be applied.': 'spearman',
            'Vector of alpha diversity values by sample.': str(aid),
            'qp-hide-metadata': 'metadata',
            'qp-hide-method': 'alpha_correlation',
            'qp-hide-paramThe correlation test to be applied.': 'method',
            'qp-hide-paramVector of alpha diversity values by '
            'sample.': 'alpha_diversity',
            'qp-hide-plugin': 'diversity'}
        self.data['command'] = dumps(
            ['qiime2', qiime2_version, 'Alpha diversity correlation'])
        self.data['parameters'] = dumps(params)

        jid = self.qclient.post(
            '/apitest/processing_job/', data=self.data)['job']
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = call_qiime2(self.qclient, jid, params, out_dir)
        self.assertEqual(msg, '')
        self.assertTrue(success)
        # only 1 element
        self.assertEqual(len(ainfo), 1)
        # and that element [0] should have this file
        exp = [(join(out_dir, 'alpha_correlation', 'visualization'), 'qzv')]
        self.assertEqual(ainfo[0].files, exp)

    def test_filter_samples(self):
        # let's test a failure
        params = {
            'If true, the samples selected by `metadata` or `where` '
            'parameters will be excluded from the filtered table instead of '
            'being retained.': False,
            'SQLite WHERE clause specifying sample metadata criteria that '
            'must be met to be included in the filtered feature table. If '
            'not provided, all samples in `metadata` that are also in the '
            'feature table will be retained.': '',
            'The feature table from which samples should be filtered.': '8',
            'The maximum number of features that a sample can have to be '
            'retained. If no value is provided this will default to infinity '
            '(i.e., no maximum feature filter will be applied).': '1',
            'The maximum total frequency that a sample can have to be '
            'retained. If no value is provided this will default to infinity '
            '(i.e., no maximum frequency filter will be applied).': '',
            'The minimum number of features that a sample must have to be '
            'retained.': '0',
            'The minimum total frequency that a sample must have to be '
            'retained.': '0',
            'qp-hide-metadata': 'metadata',
            'qp-hide-method': 'filter_samples',
            'qp-hide-paramThe feature table from which samples should be '
            'filtered.': 'table',
            'qp-hide-plugin': 'feature-table',
            'qp-hide-paramIf true, the samples selected by `metadata` or '
            '`where` parameters will be excluded from the filtered table '
            'instead of being retained.': 'exclude_ids',
            'qp-hide-paramSQLite WHERE clause specifying sample metadata '
            'criteria that must be met to be included in the filtered '
            'feature table. If not provided, all samples in `metadata` that '
            'are also in the feature table will be retained.': 'where',
            'qp-hide-paramThe maximum number of features that a sample can '
            'have to be retained. If no value is provided this will default '
            'to infinity (i.e., no maximum feature filter will be '
            'applied).': 'max_features',
            'qp-hide-paramThe maximum total frequency that a sample can have '
            'to be retained. If no value is provided this will default to '
            'infinity (i.e., no maximum frequency filter will be '
            'applied).': u'max_frequency',
            'qp-hide-paramThe minimum number of features that a sample must '
            'have to be retained.': 'min_features',
            'qp-hide-paramThe minimum total frequency that a sample must '
            'have to be retained.': 'min_frequency'}
        self.data['command'] = dumps(
            ['qiime2', qiime2_version, 'Filter samples from table'])
        self.data['parameters'] = dumps(params)

        jid = self.qclient.post(
            '/apitest/processing_job/', data=self.data)['job']
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = call_qiime2(self.qclient, jid, params, out_dir)
        self.assertEqual(msg, '')
        self.assertTrue(success)
        self.assertEqual(ainfo[0].files, [(
            join(out_dir, 'filter_samples', 'filtered_table',
                 'feature-table.biom'), 'biom')])
        self.assertEqual(ainfo[0].output_name, 'filtered_table')


if __name__ == '__main__':
    main()
