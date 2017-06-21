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
from qiita_client.util import system_call


STATE_UNIFRAC_METRICS = {
    "unweighted UniFrac": "unweighted",
    "weighted normalized UniFrac": "weighted-normalized",
    "weighted unnormalized UniFrac": "weighted-unnormalized"}


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
    artifact_id = int(parameters['i-table'])
    rarefy_level = int(parameters['p-sampling-depth'])
    artifact_info = qclient.get("/qiita_db/artifacts/%d/" % artifact_id)

    # getting just the biom file, [0] it should be only one
    to_rarefy = artifact_info['files']['biom'][0]
    qclient.update_job_step(job_id, "Step 2 of 2: Rarefying")
    b = load_table(to_rarefy)

    if not exists(out_dir):
        mkdir(out_dir)

    rarefied = b.subsample(rarefy_level)
    if rarefied.sum() == 0:
        return False, None, "Rarefaction level too high %d" % rarefy_level

    rarefied_fp = join(out_dir, 'rarefied.biom')
    with biom_open(rarefied_fp, 'w') as bf:
        rarefied.to_hdf5(bf, "Qiita's Qiime2 plugin")

    ainfo = [ArtifactInfo('o-table', 'BIOM', [(rarefied_fp, 'biom')])]

    return True, ainfo, ""


def beta_diversity(qclient, job_id, parameters, out_dir):
    """generate beta diversity calculations

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
    out_dir = join(out_dir, 'beta_diversity')
    if not exists(out_dir):
        mkdir(out_dir)

    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    artifact_id = parameters['i-table']
    metric = parameters['p-metric']
    tree = parameters['i-tree']
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    biom_fpi = artifact_info['files']['biom'][0]
    biom_qza = join(out_dir, 'q2-biom.qza')

    qclient.update_job_step(
        job_id, "Step 2 of 4: Converting Qiita artifacts to Q2 artifact")
    # converting biom
    cmd = ('qiime tools import --input-path %s --output-path %s '
           '--type "FeatureTable[Frequency]' % (biom_fpi, biom_qza))
    b = load_table(biom_fpi)
    counts = list(map(sum, b.iter_data()))
    if min(counts) == max(counts):
        cmd += " % Properties(['uniform-sampling'])\""
    else:
        cmd += '"'
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error converting biom:\nStd out: %s\nStd err: %s"
                     % (std_out, std_err))
        return False, None, error_msg
    # converting tree
    if tree is not None:
        qza_tree = join(out_dir, 'tree.qza')
        cmd = ('qiime tools import --input-path %s --type Phylogeny[Rooted] '
               '--output-path %s' % (tree, qza_tree))
        tree = qza_tree
        std_out, std_err, return_value = system_call(cmd)
        if return_value != 0:
            error_msg = ("Error converting biom:\nStd out: %s\nStd err: %s"
                         % (std_out, std_err))
            return False, None, error_msg

    qclient.update_job_step(
        job_id, "Step 3 of 4: Calculating beta diversity: %s" % (metric))
    if tree is not None and metric in STATE_UNIFRAC_METRICS:
        su_metric = STATE_UNIFRAC_METRICS[metric]
        dtx_fp = join(out_dir, '%s.qza' % su_metric)
        cmd = ('qiime state-unifrac %s --i-table %s --i-phylogeny %s '
               '--o-distance-matrix %s' % (su_metric, biom_qza, tree, dtx_fp))
    elif metric not in STATE_UNIFRAC_METRICS and tree is None:
        dtx_fp = join(out_dir, '%s.qza' % metric)
        cmd = ('qiime diversity beta --i-table %s --p-metric %s '
               '--o-distance-matrix %s' % (biom_qza, metric, dtx_fp))
    else:
        return False, None, ('Phylogentic metric %s selected but no tree '
                             'exists' % metric)

    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error in beta div %s:\nStd out: %s\nStd err: %s"
                     % (metric, std_out, std_err))
        return False, None, error_msg

    qclient.update_job_step(
        job_id, "Step 4 of 4: Converting Q2 to Qiita artifacts")
    fdir = join(out_dir, 'dtx')
    ffp = join(fdir, 'distance-matrix.tsv')
    cmd = "qiime tools export --output-dir %s %s" % (fdir, dtx_fp)
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error in Q2 -> Qiita conversion:\nStd out: "
                     "%s\nStd err: %s" % (std_out, std_err))
        return False, None, error_msg

    ainfo = [ArtifactInfo('distance-matrix', 'distance_matrix',
                          [(ffp, 'plain_text')])]
    return True, ainfo, ""
