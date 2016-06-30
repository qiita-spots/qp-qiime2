# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from tempfile import mkdtemp
from os import environ, remove
from os.path import exists, isdir
from shutil import rmtree
from json import dumps

from qiita_client import QiitaClient
from biom import load_table
import numpy as np
import numpy.testing as npt

from qp_qiime2.plugin import execute_job, execute_task

CLIENT_ID = '19ndkO3oMKsoChjVVWluF7QkxHRfYhTKSFbAVt8IhK7gZgDaO4'
CLIENT_SECRET = ('J7FfQ7CQdOxuKhQAf1eoGgBAE81Ns8Gu3EKaWFm3IO2JKh'
                 'AmmCWZuabe0O5Mp28s1')


class PluginTests(TestCase):
    @classmethod
    def setUpClass(cls):
        server_cert = environ.get('QIITA_SERVER_CERT', None)
        cls.qclient = QiitaClient("https://localhost:21174", CLIENT_ID,
                                  CLIENT_SECRET, server_cert=server_cert)

    @classmethod
    def tearDownClass(cls):
        cls.qclient.post("/apitest/reset/")

    def setUp(self):
        self.out_dir = mkdtemp()
        self._clean_up_files = [self.out_dir]

        data = {'command': 8,
                'parameters': dumps({'table': 4, 'depth': 100}),
                'status': 'queued'}
        self.job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_execute_task(self):
        job_info = self.qclient.get_job_info(self.job_id)
        success, artifacts_info, error_msg = execute_task(
            self.qclient, self.job_id, job_info)
        obs = load_table(artifacts_info[0].files[0][0]).sum(axis='sample')
        exp = np.array([100] * 7)
        npt.assert_almost_equal(obs, exp)

    def test_execute_job(self):
        execute_job("https://localhost:21174", self.job_id, self.out_dir)

        obs = self.qclient.get_job_info(self.job_id)
        self.assertEqual(obs['status'], 'success')

if __name__ == '__main__':
    main()
