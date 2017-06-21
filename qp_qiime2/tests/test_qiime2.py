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
from qp_qiime2.qiime2 import (rarefy, beta_diversity)


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
                'command': dumps(['qiime2', qiime2_version, 'Rarefy']),
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
                'command': dumps(['qiime2', qiime2_version, 'Rarefy']),
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
                'command': dumps(['qiime2', qiime2_version, 'Rarefy']),
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
            'i-table': aid, 'p-metric': 'euclidean',
            'i-tree': None}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version, 'beta_diversity']),
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
                'distance_matrix')]
        self.assertEqual(ainfo[0].files, exp)

        params['p-metric'] = 'unweighted'
        params['i-tree'] = join(
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
                'distance_matrix')]
        self.assertEqual(ainfo[0].files, exp)

        # To avoid having to set up all these files, we are gonna test
        # that if phylogentic and no tree it fails
        params['i-tree'] = None
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = beta_diversity(
            self.qclient, jid, params, out_dir)
        self.assertFalse(success)
        self.assertEqual(msg, 'Phylogentic metric unweighted selected but '
                              'no tree exists')

    def test_beta_errors(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # no rarefied - this testes that the non rarefied table conversion
        # works and that if qiime fails it raises the correct error_msg
        params = {
            'i-table': 5, 'p-metric': 'braycurtis', 'i-tree': None}
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qiime2', qiime2_version, 'beta_diversity']),
                'status': 'running',
                'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        success, ainfo, msg = beta_diversity(
            self.qclient, jid, params, out_dir)

        self.assertIn("Argument to input 'table' is not a subtype of "
                      "FeatureTable[Frequency]", msg)
        self.assertFalse(success)


if __name__ == '__main__':
    main()
