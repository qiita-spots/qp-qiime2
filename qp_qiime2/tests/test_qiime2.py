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
from os.path import exists, isdir, join
from biom import load_table

from qiita_client.testing import PluginTestCase

from qiime2 import __version__ as qiime2_version

from qp_qiime2 import plugin
from qp_qiime2.qiime2 import (rarefy)


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
        self.assertEqual(True, success)
        self.assertEqual('', msg)
        self.assertEqual([(join(out_dir, 'rarefy', 'rarefied.biom'), 'biom')],
                         ainfo[0].files)
        self.assertEqual('rarefied table @ 2', ainfo[0].output_name)

        # testing that the table is actually rarefied, [0] cause there is only
        # one element, and [0][0] from that element we want the first element
        # of the first tuple
        rb = load_table(ainfo[0].files[0][0])
        # 2 * 7 cause we rarefied at 2 sequences per sample and we have 7
        # samples
        self.assertEqual(2 * 7, rb.sum())

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
        self.assertEqual(False, success)
        self.assertIsNone(ainfo)
        self.assertEqual('Rarefaction level too high 200000', msg)


if __name__ == '__main__':
    main()
