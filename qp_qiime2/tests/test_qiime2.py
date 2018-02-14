# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from sys import maxsize
from os import remove
from shutil import rmtree
from tempfile import mkdtemp
from json import dumps
from os.path import exists, isdir, join, realpath, dirname
from biom import load_table

from qiita_client.testing import PluginTestCase

from qiime2 import __version__ as qiime2_version
from q2_diversity.plugin_setup import plugin as q2div_plugin

from qp_qiime2 import plugin
from qp_qiime2.qiime2 import (
    rarefy, beta_diversity, pcoa, beta_correlation, alpha_diversity,
    alpha_correlation, taxa_barplot, filter_samples, emperor,
    beta_group_significance, BETA_DIVERSITY_METRICS, STATE_UNIFRAC_METRICS,
    ALPHA_DIVERSITY_METRICS, ALPHA_CORRELATION_METHODS,
    BETA_CORRELATION_METHODS, BETA_GROUP_SIG_METHODS)


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
        params = {'Sampling depth': 2, 'BIOM table': 5}
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
        self.assertEqual(ainfo[0].output_name, 'Rarefied table')

        # testing that the table is actually rarefied, [0] cause there is only
        # one element, and [0][0] from that element we want the first element
        # of the first tuple
        rb = load_table(ainfo[0].files[0][0])
        # 2 * 7 cause we rarefied at 2 sequences per sample and we have 7
        # samples
        self.assertEqual(rb.sum(), 2 * 7)

    def test_rarefy_error(self):
        params = {'Sampling depth': 200000, 'BIOM table': 5}
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
        params = {'Sampling depth': 10, 'BIOM table': 5}
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
            'Phylogenetic tree': 'None', 'Number of jobs': 1}
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
        params['Adjust variance (phylogenetic only)'] = False
        params['Bypass tips (phylogenetic only)'] = False
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
        params['Diversity metric'] = 'Euclidean distance'
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = beta_diversity(
            self.qclient, jid, params, out_dir)
        self.assertFalse(success)
        self.assertEqual(msg, 'Error. Metric: euclidean (is '
                              'phylogenetic: %s), tree: None' %
                              params['Phylogenetic tree'])

        params['Phylogenetic tree'] = "None"
        params['Diversity metric'] = 'Unweighted UniFrac'
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = beta_diversity(
          self.qclient, jid, params, out_dir)
        self.assertFalse(success)
        self.assertEqual(msg, 'Error. Metric: unweighted UniFrac (is '
                              'phylogenetic: True), tree: None')

    def test_pcoa(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # qiime2 currently only works with rarefied tables so we need to
        # rarefy it
        params = {'Sampling depth': 10, 'BIOM table': 5}
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
            'Phylogenetic tree': 'None', 'Number of jobs': 1}
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
        params = {'Distance matrix': aid}
        data = {'user': 'demo@microbio.me',
                'command': dumps(
                    ['qiime2', qiime2_version,
                     'Perform Principal Coordinates Analysis (PCoA)']),
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
        params = {'Sampling depth': 10, 'BIOM table': 8}
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
            'Phylogenetic tree': 'None', 'Number of jobs': 1}
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
        params = {'Distance matrix': aid,
                  'Metadata category': 'samp_salinity',
                  'Correlation method': 'Spearman',
                  'Number of permutations': 5}
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

        # testing failure here, just to avoid reduplicating all the code above
        params = {'Distance matrix': aid,
                  'Metadata category': 'common_name',
                  'Correlation method': 'Spearman',
                  'Number of permutations': 5}
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
        params = {'Sampling depth': 10, 'BIOM table': 5}
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
        params['Diversity metric'] = 'Number of distinct features'
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = alpha_diversity(
            self.qclient, jid, params, out_dir)
        self.assertFalse(success)
        self.assertEqual(msg, 'Error. Metric: observed_otus (is '
                              'phylogenetic: False), tree: %s' %
                              params['Phylogenetic tree'])

        params['Phylogenetic tree'] = "None"
        params['Diversity metric'] = "Faith's Phylogenetic Diversity"
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = alpha_diversity(
            self.qclient, jid, params, out_dir)
        self.assertFalse(success)
        self.assertEqual(msg, 'Error. Metric: faith_pd (is phylogenetic: '
                              'True), tree: None')

    def test_alpha_correlation(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # qiime2 currently only works with rarefied tables so we need to
        # rarefy it
        params = {'Sampling depth': 10, 'BIOM table': 8}
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
        params = {'Alpha vectors': aid, 'Correlation method': 'Spearman'}
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

        params = {'BIOM table': 8}
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
            'BIOM table': 8,
            'Minimum feature frequency across samples': 1,
            'Maximum feature frequency across samples': maxsize,
            'Minimum features per sample': 1,
            'Maximum features per sample': maxsize,
            'SQLite WHERE-clause': '',
            'Exclude ids selected by where parameter': False
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
        exp = [(join(out_dir, 'filter_samples/filtered.biom'), 'biom')]
        self.assertEqual(ainfo[0].files, exp)
        self.assertIn(
            ainfo[0].output_name,
            plugin.task_dict['Filter samples by metadata'].outputs)

        b = load_table(exp[0][0])
        self.assertIn('taxonomy', b.metadata(b.ids(axis='observation')[0],
                                             axis='observation'))
        self.assertEqual(b.metadata('3862157', axis='observation')['taxonomy'],
                         ['k__Bacteria', 'p__Proteobacteria',
                          'c__Alphaproteobacteria', 'o__', 'f__',
                          'g__', 's__'])
        self.assertEqual(b.metadata('185100', axis='observation')['taxonomy'],
                         ['k__Bacteria', 'p__Proteobacteria',
                          'c__Deltaproteobacteria', 'o__Bdellovibrionales',
                          'f__Bacteriovoracaceae', 'g__', 's__'])

    def test_emperor(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # qiime2 currently only works with rarefied tables so we need to
        # rarefy it
        params = {'Sampling depth': 10, 'BIOM table': 8}
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
            'Phylogenetic tree': 'None', 'Number of jobs': 1}
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
        params = {'Distance matrix': aid}
        data = {'user': 'demo@microbio.me',
                'command': dumps(
                    ['qiime2', qiime2_version,
                     'Perform Principal Coordinates Analysis (PCoA)']),
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
        params = {'Ordination results': aid, 'Custom axis': 'latitude'}
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
        params = {'Sampling depth': 10, 'BIOM table': 8}
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
            'Phylogenetic tree': 'None', 'Number of jobs': 1}
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
        params = {'Distance matrix': aid,
                  'Metadata category': 'samp_salinity',
                  'Comparison type': 'Pairwise',
                  'Method': 'PERMANOVA',
                  'Number of permutations': 5}
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

    def test_metrics(self):
        # Test beta diversity metrics
        beta_methods = q2div_plugin.methods['beta'].signature.parameters[
            'metric'].qiime_type.predicate.choices
        beta_alt_methods = q2div_plugin.methods[
            'beta_phylogenetic_alt'].signature.parameters[
                'metric'].qiime_type.predicate.choices
        q2_metrics = beta_methods.union(beta_alt_methods)
        qp_metrics = set(BETA_DIVERSITY_METRICS.values()).union(
            STATE_UNIFRAC_METRICS.values()).difference(STATE_UNIFRAC_METRICS)
        self.assertEqual(q2_metrics, qp_metrics)

        # Test alpha diversity metrics
        alpha_methods = q2div_plugin.methods['alpha'].signature.parameters[
            'metric'].qiime_type.predicate.choices
        alpha_phyl_methods = q2div_plugin.methods[
            'alpha_phylogenetic'].signature.parameters[
                'metric'].qiime_type.predicate.choices
        q2_metrics = alpha_methods.union(alpha_phyl_methods)
        qp_metrics = set(ALPHA_DIVERSITY_METRICS.values())
        self.assertEqual(q2_metrics, qp_metrics)

        # Alpha correlation methods
        q2_methods = q2div_plugin.visualizers[
            'alpha_correlation'].signature.parameters[
                'method'].qiime_type.predicate.choices
        qp_methods = set(ALPHA_CORRELATION_METHODS.values())
        self.assertEqual(q2_methods, qp_methods)

        # Beta correlation methods
        q2_methods = q2div_plugin.visualizers[
            'mantel'].signature.parameters[
                'method'].qiime_type.predicate.choices
        qp_methods = set(BETA_CORRELATION_METHODS.values())
        self.assertEqual(q2_methods, qp_methods)

        # Beta group significance methods
        q2_methods = q2div_plugin.visualizers[
            'beta_group_significance'].signature.parameters[
                'method'].qiime_type.predicate.choices
        qp_methods = set(BETA_GROUP_SIG_METHODS.values())
        self.assertEqual(q2_methods, qp_methods)


if __name__ == '__main__':
    main()
