# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os import mkdir
from os.path import join, exists

import pandas as pd

from biom import load_table
from biom.util import biom_open
from qiita_client import ArtifactInfo
from qiita_client.util import system_call


STATE_UNIFRAC_METRICS = {
    "unweighted UniFrac": "unweighted",
    "weighted normalized UniFrac": "weighted-normalized",
    "weighted unnormalized UniFrac": "weighted-unnormalized"}

ALPHA_PHYLOGENETIC_METRICS = {'faith_pd'}


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

    biom_qza, metric, tree = _diversity_init_steps(
        qclient, job_id, parameters, out_dir)

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

    ainfo = [ArtifactInfo('distance_matrix', 'distance_matrix',
                          [(ffp, 'plain_text')])]
    return True, ainfo, ""


def pcoa(qclient, job_id, parameters, out_dir):
    """generate pcoa calculations

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
    out_dir = join(out_dir, 'pcoa')
    if not exists(out_dir):
        mkdir(out_dir)

    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    artifact_id = parameters['i-distance-matrix']
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    dm_fp = artifact_info['files']['plain_text'][0]
    dm_qza = join(out_dir, 'q2-distance.qza')
    pcoa_qza = join(out_dir, 'q2-pcoa.qza')

    qclient.update_job_step(
        job_id, "Step 2 of 4: Converting Qiita artifacts to Q2 artifact")
    cmd = ('qiime tools import --input-path %s --output-path %s '
           '--type "DistanceMatrix"' % (dm_fp, dm_qza))
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error converting distance matrix:\nStd out: %s\n"
                     "Std err: %s" % (std_out, std_err))
        return False, None, error_msg

    qclient.update_job_step(
        job_id, "Step 3 of 4: Calculating pcoa")
    cmd = ('qiime diversity pcoa --i-distance-matrix %s --o-pcoa %s' % (
        dm_qza, pcoa_qza))

    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error in PCoA\nStd out: %s\nStd err: %s"
                     % (std_out, std_err))
        return False, None, error_msg

    qclient.update_job_step(
        job_id, "Step 4 of 4: Converting Q2 to Qiita artifacts")
    fdir = join(out_dir, 'pcoa')
    ffp = join(fdir, 'ordination.txt')
    cmd = "qiime tools export --output-dir %s %s" % (fdir, pcoa_qza)
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error in Q2 -> Qiita conversion:\nStd out: "
                     "%s\nStd err: %s" % (std_out, std_err))
        return False, None, error_msg

    ainfo = [ArtifactInfo('o-pcoa', 'ordination_results',
                          [(ffp, 'plain_text')])]
    return True, ainfo, ""


def beta_correlation(qclient, job_id, parameters, out_dir):
    """generate beta correlation calculations

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
    out_dir = join(out_dir, 'beta_correlation')
    if not exists(out_dir):
        mkdir(out_dir)

    qclient.update_job_step(job_id, "Step 1 of 3: Collecting information")
    artifact_id = parameters['i-distance-matrix']
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    dm_fp = artifact_info['files']['plain_text'][0]
    dm_qza = join(out_dir, 'q2-distance.qza')
    analysis_id = artifact_info['analysis']
    metadata = qclient.get(
        "/qiita_db/analysis/%s/metadata/" % str(analysis_id))
    metadata = pd.DataFrame.from_dict(metadata, orient='index')
    metadata_fp = join(out_dir, 'metadata.txt')
    metadata.to_csv(metadata_fp, sep='\t')
    m_metadata_category = parameters['m-metadata-category']
    p_method = parameters['p-method']
    p_permutations = parameters['p-permutations']
    o_visualization = join(out_dir, 'beta_correlation.qzv')

    qclient.update_job_step(
        job_id, "Step 2 of 3: Converting Qiita artifacts to Q2 artifact")
    cmd = ('qiime tools import --input-path %s --output-path %s '
           '--type "DistanceMatrix"' % (dm_fp, dm_qza))
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error converting distance matrix:\nStd out: %s\n"
                     "Std err: %s" % (std_out, std_err))
        return False, None, error_msg

    qclient.update_job_step(
        job_id, "Step 3 of 3: Calculating beta correlation")
    cmd = ('qiime diversity beta-correlation --i-distance-matrix %s '
           '--m-metadata-file %s --m-metadata-category %s --p-method %s '
           '--p-permutations %s --o-visualization %s' % (
               dm_qza, metadata_fp, m_metadata_category, p_method,
               p_permutations, o_visualization))

    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error in Beta Correlation\nStd out: %s\nStd err: %s"
                     % (std_out, std_err))
        return False, None, error_msg

    ainfo = [ArtifactInfo('o-visualization', 'q2_visualization',
                          [(o_visualization, 'qiime2-visualization')])]
    return True, ainfo, ""


def alpha_diversity(qclient, job_id, parameters, out_dir):
    """generate alpha diversity calculations

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
    out_dir = join(out_dir, 'alpha_diversity')
    if not exists(out_dir):
        mkdir(out_dir)

    biom_qza, metric, tree = _diversity_init_steps(
        qclient, job_id, parameters, out_dir)

    qclient.update_job_step(
        job_id, "Step 3 of 4: Calculating alpha diversity: %s" % (metric))
    alpha_fp = join(out_dir, '%s.qza' % metric)
    if tree is not None and metric in ALPHA_PHYLOGENETIC_METRICS:
        cmd = 'qiime diversity alpha-phylogenetic --i-phylogeny %s ' % tree
    elif metric not in ALPHA_PHYLOGENETIC_METRICS and tree is None:
        cmd = 'qiime diversity alpha '
    else:
        return False, None, ('Phylogentic metric %s selected but no tree '
                             'exists' % metric)
    cmd += '--i-table %s --p-metric %s --o-alpha-diversity %s' % (
        biom_qza, metric, alpha_fp)

    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error in alpha div %s:\nStd out: %s\nStd err: %s"
                     % (metric, std_out, std_err))
        return False, None, error_msg

    qclient.update_job_step(
        job_id, "Step 4 of 4: Converting Q2 to Qiita artifacts")
    fdir = join(out_dir, 'alpha')
    ffp = join(fdir, 'alpha-diversity.tsv')
    cmd = "qiime tools export --output-dir %s %s" % (fdir, alpha_fp)
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error in Q2 -> Qiita conversion:\nStd out: "
                     "%s\nStd err: %s" % (std_out, std_err))
        return False, None, error_msg

    ainfo = [ArtifactInfo('o-alpha-diversity', 'alpha_vector',
                          [(ffp, 'plain_text')])]
    return True, ainfo, ""


def _diversity_init_steps(qclient, job_id, parameters, out_dir):
    """helper function to avoid duplication of code"""

    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    artifact_id = parameters['i-table']
    metric = parameters['p-metric']
    tree = parameters['i-tree']
    if tree == 'None':
        tree = None
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

    return biom_qza, metric, tree


def alpha_correlation(qclient, job_id, parameters, out_dir):
    """generate alpha correlation calculations

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
    out_dir = join(out_dir, 'alpha_correlation')
    if not exists(out_dir):
        mkdir(out_dir)

    qclient.update_job_step(job_id, "Step 1 of 3: Collecting information")
    artifact_id = parameters['i-alpha-diversity']
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    dm_fp = artifact_info['files']['plain_text'][0]
    dm_qza = join(out_dir, 'q2-alpha-diversity.qza')
    analysis_id = artifact_info['analysis']
    metadata = qclient.get(
        "/qiita_db/analysis/%s/metadata/" % str(analysis_id))
    metadata = pd.DataFrame.from_dict(metadata, orient='index')
    metadata_fp = join(out_dir, 'metadata.txt')
    metadata.to_csv(metadata_fp, sep='\t')
    p_method = parameters['p-method']
    o_visualization = join(out_dir, 'alpha_correlation.qzv')

    qclient.update_job_step(
        job_id, "Step 2 of 3: Converting Qiita artifacts to Q2 artifact")
    cmd = ('qiime tools import --input-path %s --output-path %s '
           '--type "SampleData[AlphaDiversity]"' % (dm_fp, dm_qza))
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error converting distance matrix:\nStd out: %s\n"
                     "Std err: %s" % (std_out, std_err))
        return False, None, error_msg

    qclient.update_job_step(
        job_id, "Step 3 of 3: Calculating alpha correlation")
    cmd = ('qiime diversity alpha-correlation --i-alpha-diversity %s '
           '--m-metadata-file %s --p-method %s --o-visualization %s' % (
               dm_qza, metadata_fp, p_method, o_visualization))

    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error in Alpha Correlation\nStd out: %s\nStd err: %s"
                     % (std_out, std_err))
        return False, None, error_msg

    ainfo = [ArtifactInfo('o-visualization', 'q2_visualization',
                          [(o_visualization, 'qiime2-visualization')])]
    return True, ainfo, ""
