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
from biom import load_table

from qiita_client.testing import PluginTestCase

from qiime2 import __version__ as qiime2_version

from qp_qiime2 import plugin
from qp_qiime2.qiime2 import (rarefy, beta_diversity, pcoa, beta_correlation,
                              alpha_diversity, alpha_correlation, taxa_barplot,
                              filter_samples, emperor, beta_group_significance)


class qiime2Tests(PluginTestCase):
    def setUp(self):
        # this will allow us to see the full errors
        self.maxDiff = None

        plugin("https://localhost:21174", 'register', 'ignored')
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_rarefy(self):
        params = {'p-sampling-depth': 2, 'i-table': 5}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Rarefy features']),
                'status': 'running',
                'parameters': dumps(params)}

        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = rarefy(self.qclient, jid, params, out_dir)
        self.assertTrue(success)
        self.assertEqual(msg, '')
        self.assertEqual(ainfo[0].files,
                         [(join(out_dir, 'rarefy', 'rarefied.biom'), 'biom')])
        self.assertEqual(ainfo[0].output_name, 'o-table')

        # testing that the table is actually rarefied, [0] cause there is only
        # one element, and [0][0] from that element we want the first element
        # of the first tuple
        rb = load_table(ainfo[0].files[0][0])
        # 2 * 7 cause we rarefied at 2 sequences per sample and we have 7
        # samples
        self.assertEqual(rb.sum(), 2 * 7)

    def test_rarefy_error(self):
        params = {'p-sampling-depth': 200000, 'i-table': 5}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Rarefy features']),
                'status': 'running',
                'parameters': dumps(params)}

        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = rarefy(self.qclient, jid, params, out_dir)
        self.assertFalse(success)
        self.assertIsNone(ainfo)
        self.assertEqual(msg, 'Rarefaction level too high 200000')

    def test_beta(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # qiime2 currently only works with rarefied tables so we need to
        # rarefy it
        params = {'p-sampling-depth': 10, 'i-table': 5}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Rarefy features']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = rarefy(self.qclient, jid, params, out_dir)
        data = {'filepaths': dumps(ainfo[0].files), 'type': "BIOM",
                'name': "Rarefied biom", 'analysis': 1, 'data_type': '16S'}
        reply = self.qclient.post('/apitest/artifact/', data=data)
        aid = reply['artifact']

        # actually test non phylogenetic beta diversity
        params = {
            'BIOM table': aid, 'Diversity metric': 'Euclidean distance',
            'Phylogenetic tree': 'None'}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Calculate beta diversity']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = beta_diversity(
            self.qclient, jid, params, out_dir)
        self.assertEqual(msg, '')
        self.assertTrue(success)
        # only 1 element
        self.assertEqual(len(ainfo), 1)
        # and that element [0] should have this file
        exp = [(join(out_dir, 'beta_diversity/dtx/distance-matrix.tsv'),
                'plain_text')]
        self.assertEqual(ainfo[0].files, exp)

        params['Diversity metric'] = 'Unweighted UniFrac'
        params['Phylogenetic tree'] = join(
            dirname(realpath(__file__)), 'prune_97_gg_13_8.tre')
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = beta_diversity(
            self.qclient, jid, params, out_dir)
        self.assertEqual(msg, '')
        self.assertTrue(success)
        # only 1 element
        self.assertEqual(len(ainfo), 1)
        # and that element [0] should have this file
        exp = [(join(out_dir, 'beta_diversity/dtx/distance-matrix.tsv'),
                'plain_text')]
        self.assertEqual(ainfo[0].files, exp)

        # To avoid having to set up all these files, we are gonna test
        # that if phylogenetic and no tree it fails
        params['Phylogenetic tree'] = None
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = beta_diversity(
            self.qclient, jid, params, out_dir)
        self.assertFalse(success)
        self.assertEqual(msg, 'Phylogenetic metric unweighted UniFrac '
                              'selected but no tree exists')

    def test_pcoa(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # qiime2 currently only works with rarefied tables so we need to
        # rarefy it
        params = {'p-sampling-depth': 10, 'i-table': 5}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Rarefy features']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = rarefy(self.qclient, jid, params, out_dir)
        data = {'filepaths': dumps(ainfo[0].files), 'type': "BIOM",
                'name': "Rarefied biom", 'analysis': 1, 'data_type': '16S'}
        reply = self.qclient.post('/apitest/artifact/', data=data)
        aid = reply['artifact']

        # non phylogenetic beta diversity
        params = {
            'BIOM table': aid, 'Diversity metric': 'Euclidean distance',
            'Phylogenetic tree': 'None'}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Calculate beta diversity']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = beta_diversity(
            self.qclient, jid, params, out_dir)
        data = {'filepaths': dumps(ainfo[0].files), 'type': "distance_matrix",
                'name': "Non phylogenetic distance matrix", 'analysis': 1,
                'data_type': '16S'}
        reply = self.qclient.post('/apitest/artifact/', data=data)
        aid = reply['artifact']

        # pcoa
        params = {'i-distance-matrix': aid}
        data = {'user': 'demo@microbio.me',
                'command': dumps(
                    ['qiime2', qiime2_version,
                     'Generate principal coordinates analysis (PCoA)']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = pcoa(self.qclient, jid, params, out_dir)

        self.assertEqual(msg, '')
        self.assertTrue(success)
        # only 1 element
        self.assertEqual(len(ainfo), 1)
        # and that element [0] should have this file
        exp = [(join(out_dir, 'pcoa/pcoa/ordination.txt'), 'plain_text')]
        self.assertEqual(ainfo[0].files, exp)

    def test_beta_correlation(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # qiime2 currently only works with rarefied tables so we need to
        # rarefy it
        params = {'p-sampling-depth': 10, 'i-table': 8}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Rarefy features']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = rarefy(self.qclient, jid, params, out_dir)
        data = {'filepaths': dumps(ainfo[0].files), 'type': "BIOM",
                'name': "Rarefied biom", 'analysis': 1, 'data_type': '16S'}
        reply = self.qclient.post('/apitest/artifact/', data=data)
        aid = reply['artifact']

        # non phylogenetic beta diversity
        params = {
            'BIOM table': aid, 'Diversity metric': 'Euclidean distance',
            'Phylogenetic tree': 'None'}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Calculate beta diversity']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = beta_diversity(
            self.qclient, jid, params, out_dir)
        data = {'filepaths': dumps(ainfo[0].files), 'type': "distance_matrix",
                'name': "Non phylogenetic distance matrix", 'analysis': 1,
                'data_type': '16S'}
        reply = self.qclient.post('/apitest/artifact/', data=data)
        aid = reply['artifact']

        # beta_correlation
        # 1 using that analysis
        params = {'i-distance-matrix': aid,
                  'm-metadata-category': 'samp_salinity',
                  'p-method': 'spearman', 'p-permutations': 5}
        data = {'user': 'demo@microbio.me',
                'command': dumps([
                    'qiime2', qiime2_version, 'Calculate beta correlation']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = beta_correlation(
            self.qclient, jid, params, out_dir)

        self.assertEqual(msg, '')
        self.assertTrue(success)
        # only 1 element
        self.assertEqual(len(ainfo), 1)
        # and that element [0] should have this file
        exp = [(join(out_dir, 'beta_correlation/beta_correlation.qzv'), 'qzv')]
        self.assertEqual(ainfo[0].files, exp)

        # testing faillure here, just to avoid reduplicating all the code above
        params = {'i-distance-matrix': aid,
                  'm-metadata-category': 'common_name', 'p-method': 'spearman',
                  'p-permutations': 5}
        data = {'user': 'demo@microbio.me',
                'command': dumps([
                    'qiime2', qiime2_version, 'Calculate beta correlation']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = beta_correlation(
            self.qclient, jid, params, out_dir)

        self.assertIn("Unable to parse string", msg)
        self.assertFalse(success)

    def test_alpha(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # qiime2 currently only works with rarefied tables so we need to
        # rarefy it
        params = {'p-sampling-depth': 10, 'i-table': 5}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Rarefy features']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = rarefy(self.qclient, jid, params, out_dir)
        data = {'filepaths': dumps(ainfo[0].files), 'type': "BIOM",
                'name': "Rarefied biom", 'analysis': 1, 'data_type': '16S'}
        reply = self.qclient.post('/apitest/artifact/', data=data)
        aid = reply['artifact']

        # actually test non phylogenetic alpha diversity
        params = {
            'BIOM table': aid,
            'Diversity metric': 'Number of distinct features',
            'Phylogenetic tree': 'None'}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Calculate alpha diversity']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = alpha_diversity(
            self.qclient, jid, params, out_dir)
        self.assertEqual(msg, '')
        self.assertTrue(success)
        # only 1 element
        self.assertEqual(len(ainfo), 1)
        # and that element [0] should have this file
        exp = [(join(out_dir, 'alpha_diversity/alpha/alpha-diversity.tsv'),
                'plain_text')]
        self.assertEqual(ainfo[0].files, exp)

        params['Diversity metric'] = "Faith's Phylogenetic Diversity"
        params['Phylogenetic tree'] = join(
            dirname(realpath(__file__)), 'prune_97_gg_13_8.tre')
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = alpha_diversity(
            self.qclient, jid, params, out_dir)
        self.assertEqual(msg, '')
        self.assertTrue(success)
        # only 1 element
        self.assertEqual(len(ainfo), 1)
        # and that element [0] should have this file
        exp = [(join(out_dir, 'alpha_diversity/alpha/alpha-diversity.tsv'),
                'plain_text')]
        self.assertEqual(ainfo[0].files, exp)

        # To avoid having to set up all these files, we are gonna test
        # that if phylogenetic and no tree it fails
        params['Phylogenetic tree'] = None
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = alpha_diversity(
            self.qclient, jid, params, out_dir)
        self.assertFalse(success)
        self.assertEqual(msg, 'Phylogenetic metric faith_pd selected '
                              'but no tree exists')

    def test_alpha_correlation(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # qiime2 currently only works with rarefied tables so we need to
        # rarefy it
        params = {'p-sampling-depth': 10, 'i-table': 8}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Rarefy features']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = rarefy(self.qclient, jid, params, out_dir)
        data = {'filepaths': dumps(ainfo[0].files), 'type': "BIOM",
                'name': "Rarefied biom", 'analysis': 1, 'data_type': '16S'}
        reply = self.qclient.post('/apitest/artifact/', data=data)
        aid = reply['artifact']

        # non phylogenetic alpha diversity
        params = {
            'BIOM table': aid,
            'Diversity metric': 'Number of distinct features',
            'Phylogenetic tree': 'None'}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Calculate alpha diversity']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = alpha_diversity(
            self.qclient, jid, params, out_dir)
        data = {'filepaths': dumps(ainfo[0].files), 'type': "alpha_vector",
                'name': "Non phylogenetic alpha vector", 'analysis': 1,
                'data_type': '16S'}
        reply = self.qclient.post('/apitest/artifact/', data=data)
        aid = reply['artifact']

        # alpha_correlation
        # 1 using that analysis
        params = {'i-alpha-diversity': aid, 'p-method': 'spearman'}
        data = {'user': 'demo@microbio.me',
                'command': dumps([
                    'qiime2', qiime2_version, 'Calculate alpha correlation']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = alpha_correlation(
            self.qclient, jid, params, out_dir)

        self.assertEqual(msg, '')
        self.assertTrue(success)
        # only 1 element
        self.assertEqual(len(ainfo), 1)
        # and that element [0] should have this file
        exp = [(join(out_dir, 'alpha_correlation/alpha_correlation.qzv'),
               'qzv')]
        self.assertEqual(ainfo[0].files, exp)

    def test_taxa_barplot(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        params = {'i-table': 8}
        data = {'user': 'demo@microbio.me',
                'command': dumps([
                    'qiime2', qiime2_version, 'Summarize taxa']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = taxa_barplot(self.qclient, jid, params, out_dir)

        self.assertEqual(msg, '')
        self.assertTrue(success)
        # only 1 element
        self.assertEqual(len(ainfo), 1)
        # and that element [0] should have this file
        exp = [(join(out_dir, 'taxa_barplot/taxa-barplot.qzv'), 'qzv')]
        self.assertEqual(ainfo[0].files, exp)

    def test_filter_samples(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        params = {
            'i-table': '8',
            'p-min-frequency': '5',
            'p-max-frequency': '10',
            'p-min-features': '5',
            'p-max-features': '9223372036854775807',
            'p-where': ''
        }
        data = {'user': 'demo@microbio.me',
                'command': dumps([
                    'qiime2', qiime2_version, 'Filter samples by metadata']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = filter_samples(
            self.qclient, jid, params, out_dir)

        self.assertEqual(msg, '')
        self.assertTrue(success)
        # only 1 element
        self.assertEqual(len(ainfo), 1)
        # and that element [0] should have this file
        exp = [(join(out_dir,
                     'filter_samples/filter_samples/feature-table.biom'),
                'biom')]
        self.assertEqual(ainfo[0].files, exp)

    def test_emperor(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # qiime2 currently only works with rarefied tables so we need to
        # rarefy it
        params = {'p-sampling-depth': 10, 'i-table': 8}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Rarefy features']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = rarefy(self.qclient, jid, params, out_dir)
        data = {'filepaths': dumps(ainfo[0].files), 'type': "BIOM",
                'name': "Rarefied biom", 'analysis': 1, 'data_type': '16S'}
        reply = self.qclient.post('/apitest/artifact/', data=data)
        aid = reply['artifact']

        # non phylogenetic beta diversity
        params = {
            'i-table': aid, 'p-metric': 'euclidean',
            'i-tree': 'None'}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Calculate beta diversity']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = beta_diversity(
            self.qclient, jid, params, out_dir)
        data = {'filepaths': dumps(ainfo[0].files), 'type': "distance_matrix",
                'name': "Non phylogenetic distance matrix", 'analysis': 1,
                'data_type': '16S'}
        reply = self.qclient.post('/apitest/artifact/', data=data)
        aid = reply['artifact']

        # pcoa
        params = {'i-distance-matrix': aid}
        data = {'user': 'demo@microbio.me',
                'command': dumps(
                    ['qiime2', qiime2_version,
                     'Generate principal coordinates analysis (PCoA)']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = pcoa(self.qclient, jid, params, out_dir)
        data = {'filepaths': dumps(ainfo[0].files),
                'type': "ordination_results", 'name': "Non phylogenetic PCoA",
                'analysis': 1, 'data_type': '16S'}
        reply = self.qclient.post('/apitest/artifact/', data=data)
        aid = reply['artifact']

        # emperor
        # 1 using that analysis
        params = {'i-pcoa': aid, 'p-custom-axis': 'latitude'}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Custom-axis Emperor plot']),
                'status': 'running', 'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = emperor(self.qclient, jid, params, out_dir)

        self.assertEqual(msg, '')
        self.assertTrue(success)
        # only 1 element
        self.assertEqual(len(ainfo), 1)
        # and that element [0] should have this file
        exp = [(join(out_dir, 'emperor/q2-emperor.qzv'), 'qzv')]
        self.assertEqual(ainfo[0].files, exp)

    def test_beta_group_significance(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # qiime2 currently only works with rarefied tables so we need to
        # rarefy it
        params = {'p-sampling-depth': 10, 'i-table': 8}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Rarefy features']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = rarefy(self.qclient, jid, params, out_dir)
        data = {'filepaths': dumps(ainfo[0].files), 'type': "BIOM",
                'name': "Rarefied biom", 'analysis': 1, 'data_type': '16S'}
        reply = self.qclient.post('/apitest/artifact/', data=data)
        aid = reply['artifact']

        # non phylogenetic beta diversity
        params = {
            'i-table': aid, 'p-metric': 'euclidean',
            'i-tree': 'None'}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version,
                                  'Calculate beta diversity']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = beta_diversity(
            self.qclient, jid, params, out_dir)
        data = {'filepaths': dumps(ainfo[0].files), 'type': "distance_matrix",
                'name': "Non phylogenetic distance matrix", 'analysis': 1,
                'data_type': '16S'}
        reply = self.qclient.post('/apitest/artifact/', data=data)
        aid = reply['artifact']

        # beta_group_significance
        # 1 using that analysis
        params = {'i-distance-matrix': aid,
                  'm-metadata-category': 'samp_salinity',
                  'p-pairwise': 'p-pairwise',
                  'p-method': 'permanova', 'p-permutations': 5}
        data = {'user': 'demo@microbio.me',
                'command': dumps([
                    'qiime2', qiime2_version,
                    'Calculate beta group significance']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = beta_group_significance(
            self.qclient, jid, params, out_dir)

        self.assertEqual(msg, '')
        self.assertTrue(success)
        # only 1 element
        self.assertEqual(len(ainfo), 1)
        # and that element [0] should have this file
        exp = [(join(out_dir,
                     'beta_group_significance/beta_group_significance.qzv'),
                'qzv')]
        self.assertEqual(ainfo[0].files, exp)


if __name__ == '__main__':
    main()
