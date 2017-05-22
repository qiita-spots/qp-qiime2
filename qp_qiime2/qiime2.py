# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os import mkdir
from os.path import join, exists

from biom import load_table
from biom.util import biom_open
from qiita_client import ArtifactInfo


def rarefy(qclient, job_id, parameters, out_dir):
    """rarefy a table

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to rarefy
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    boolean, list, str
        The results of the job
    """
    out_dir = join(out_dir, 'rarefy')

    qclient.update_job_step(job_id, "Step 1 of 2: Collecting information")
    artifact_id = parameters['i-table']
    rarefy_level = parameters['p-sampling-depth']
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)

    # getting just the biom file, [0] it should be only one
    to_rarefy = artifact_info['files']['biom'][0]
    qclient.update_job_step(job_id, "Step 2 of 2: Rarefying")
    b = load_table(to_rarefy)

    if not exists(out_dir):
        mkdir(out_dir)

    rarefied = b.subsample(rarefy_level)
    rarefied_fp = join(out_dir, 'rarefied.biom')
    with biom_open(rarefied_fp, 'w') as bf:
        rarefied.to_hdf5(bf, "Qiita's Qiime2 plugin")

    ainfo = [ArtifactInfo('rarefied table @ %d' % rarefy_level, 'BIOM',
                          [(rarefied_fp, 'biom')])]

    return True, ainfo, ""
