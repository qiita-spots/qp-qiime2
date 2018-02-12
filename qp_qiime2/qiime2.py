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
    "unweighted UniFrac": "unweighted_unifrac",
    "weighted normalized UniFrac": "weighted_normalized_unifrac",
    "weighted unnormalized UniFrac": "weighted_unifrac",
    "generalized UniFrac": "generalized_unifrac"}

ALPHA_PHYLOGENETIC_METRICS = {'faith_pd'}

ALPHA_DIVERSITY_METRICS = {
    "Abundance-based Coverage Estimator (ACE) metric": "ace",
    "Berger-Parker dominance index": "berger_parker_d",
    "Brillouin's index": "brillouin_d",
    "Chao1 index": "chao1",
    "Chao1 confidence interval": "chao1_ci",
    "Dominance measure": "dominance",
    "Effective number of species (ENS)/"
    "Probability of intra-or interspecific encounter (PIE) metric": "enspie",
    "Esty's confidence interval": "esty_ci",
    "Faith's Phylogenetic Diversity": "faith_pd",
    "Fisher's index": "fisher_alpha",
    "Gini index": "gini_index",
    "Good's coverage of counts": "goods_coverage",
    "Heip's evenness measure": "heip_e",
    "Kempton-Taylor Q index": "kempton_taylor_q",
    "Lladser's confidence interval": "lladser_ci",
    "Lladser's point estimate": "lladser_pe",
    "Margalef's richness index": "margalef",
    "McIntosh dominance index D": "mcintosh_d",
    "McIntosh evenness index E": "mcintosh_e",
    "Menhinick's richness index": "menhinick",
    "Michaelis-Menten fit to rarefaction curve of obeserved OTUs":
        "michaelis_menten_fit",
    "Number of distinct features": "observed_otus",
    "Number of double occurrences": "doubles",
    "Number of single occurrences": "singles",
    "Number of observed features, including singles and doubles": "osd",
    "Pielou's evenness": "pielou_e",
    "Robbins' estimator": "robbins",
    "Shannon's index": "shannon",
    "Simpson's index": "simpson",
    "Simpson's evenness measure E": "simpson_e",
    "Strong's dominance index (Dw)": "strong"}

ALPHA_CORRELATION_METHODS = {
    "Spearman": "spearman",
    "Pearson": "pearson"}


BETA_DIVERSITY_METRICS = {
    "Bray-Curtis dissimilarity": "braycurtis",
    "Canberra distance": "canberra",
    "Chebysev distance": "chebyshev",
    "City-block distance": "cityblock",
    "Correlation coefficient": "correlation",
    "Cosine similarity": "cosine",
    "Dice measure": "dice",
    "Euclidean distance": "euclidean",
    "Hamming distance": "hamming",
    "Jaccard similarity index": "jaccard",
    "Kulczynski dissimilarity index": "kulsinski",
    "Mahalanobis distance": "mahalanobis",
    "Matching components": "matching",
    "Rogers-Tanimoto distance": "rogerstanimoto",
    "Russell-Rao coefficients": "russellrao",
    "Sokal-Michener coefficient": "sokalmichener",
    "Sokal-Sneath index": "sokalsneath",
    "Species-by-species Euclidean": "seuclidean",
    "Squared Euclidean": "sqeuclidean",
    "Weighted Minkowski metric": "wminkowski",
    "Yule index": "yule",
    "Unweighted UniFrac": "unweighted UniFrac",
    "Weighted normalized UniFrac": "weighted normalized UniFrac",
    "Weighted unnormalized UniFrac": "weighted unnormalized UniFrac",
    "Generalized UniFrac": "generalized UniFrac"}

BETA_CORRELATION_METHODS = {
    "Spearman": "spearman",
    "Pearson": "pearson"}

BETA_GROUP_SIG_METHODS = {
    "PERMANOVA": "permanova",
    "ANOSIM": "anosim"}

BETA_GROUP_SIG_TYPE = {
    "Pairwise": "p-pairwise",
    "Non-pairwise": "p-no-pairwise"}


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
    artifact_id = int(parameters['BIOM table'])
    rarefy_level = int(parameters['Sampling depth'])
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

    ainfo = [ArtifactInfo('Rarefied table', 'BIOM', [(rarefied_fp, 'biom')])]

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
        The parameter values for beta diversity
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
    artifact_id = parameters['BIOM table']
    metric = BETA_DIVERSITY_METRICS[parameters['Diversity metric']]
    tree = parameters['Phylogenetic tree']
    if tree == 'None':
        tree = None
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    biom_fpi = artifact_info['files']['biom'][0]
    biom_qza = join(out_dir, 'q2-biom.qza')
    num_jobs = parameters['Number of jobs']

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
            error_msg = ("Error converting tree:\nStd out: %s\nStd err: %s"
                         % (std_out, std_err))
            return False, None, error_msg

    qclient.update_job_step(
        job_id, "Step 3 of 4: Calculating beta diversity: %s" % (metric))
    if tree is not None and metric in STATE_UNIFRAC_METRICS:
        su_metric = STATE_UNIFRAC_METRICS[metric]
        dtx_fp = join(out_dir, '%s.qza' % su_metric)
        cmd = ('qiime diversity beta-phylogenetic-alt --p-metric %s '
               '--i-table %s --i-phylogeny %s --o-distance-matrix %s '
               '--p-n-jobs %s'
               % (su_metric, biom_qza, tree, dtx_fp, num_jobs))
        if parameters['Adjust variance (phylogenetic only)']:
            cmd += ' --p-variance-adjusted'
        if parameters['Bypass tips (phylogenetic only)']:
            cmd += ' --p-bypass-tips'
        if su_metric == 'generalized_unifrac':
            cmd += '--p-alpha %s' % parameters[
                'Alpha value (Generalized Unifrac only)']
    elif metric not in STATE_UNIFRAC_METRICS and tree is None:
        dtx_fp = join(out_dir, '%s.qza' % metric)
        cmd = ('qiime diversity beta --i-table %s --p-metric %s '
               '--o-distance-matrix %s --p-n-jobs %s'
               % (biom_qza, metric, dtx_fp, num_jobs))
    else:
        return False, None, ('Phylogenetic metric %s selected but no tree '
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

    ainfo = [ArtifactInfo('Distance matrix', 'distance_matrix',
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
        The parameter values for pcoa
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
    artifact_id = parameters['Distance matrix']
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

    ainfo = [ArtifactInfo('Ordination results', 'ordination_results',
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
        The parameter values for beta correlation
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

    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    artifact_id = parameters['Distance matrix']
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    dm_fp = artifact_info['files']['plain_text'][0]
    dm_qza = join(out_dir, 'q2-distance.qza')
    metadata_qza = join(out_dir, 'q2-metadata.qza')
    analysis_id = artifact_info['analysis']
    metadata = qclient.get(
        "/qiita_db/analysis/%s/metadata/" % str(analysis_id))
    metadata = pd.DataFrame.from_dict(metadata, orient='index')
    metadata_fp = join(out_dir, 'metadata.txt')
    metadata.to_csv(metadata_fp, sep='\t')
    m_metadata_category = parameters['Metadata category']
    p_method = BETA_CORRELATION_METHODS[parameters['Correlation method']]
    p_permutations = parameters['Number of permutations']
    o_visualization = join(out_dir, 'beta_correlation.qzv')

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
        job_id, "Step 3 of 4: Calculating distance matrix from metadata "
        "category")
    cmd = ('qiime metadata distance-matrix --m-metadata-file %s '
           '--m-metadata-category %s --o-distance-matrix %s' % (
               metadata_fp, m_metadata_category, metadata_qza))
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error calculating distance matrix from metadata:\n"
                     "Std out: %s\nStd err: %s" % (std_out, std_err))
        return False, None, error_msg

    qclient.update_job_step(
        job_id, "Step 4 of 4: Calculating beta correlation")
    cmd = ('qiime diversity mantel --i-dm1 %s --i-dm2 %s --p-method %s '
           '--p-permutations %s --p-intersect-ids '
           '--p-label1 "Distance Matrix" --p-label2 "%s" '
           '--o-visualization %s' % (
              dm_qza, metadata_qza, p_method, p_permutations,
              m_metadata_category, o_visualization))

    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error in Beta Correlation\nStd out: %s\nStd err: %s"
                     % (std_out, std_err))
        return False, None, error_msg

    ainfo = [ArtifactInfo('Beta correlation visualization', 'q2_visualization',
                          [(o_visualization, 'qzv')])]
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
        The parameter values for alpha diversity
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

    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    artifact_id = parameters['BIOM table']
    metric = ALPHA_DIVERSITY_METRICS[parameters['Diversity metric']]
    tree = parameters['Phylogenetic tree']
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
            error_msg = ("Error converting tree:\nStd out: %s\nStd err: %s"
                         % (std_out, std_err))
            return False, None, error_msg

    qclient.update_job_step(
        job_id, "Step 3 of 4: Calculating alpha diversity: %s" % (metric))
    alpha_fp = join(out_dir, '%s.qza' % metric)
    if tree is not None and metric in ALPHA_PHYLOGENETIC_METRICS:
        cmd = 'qiime diversity alpha-phylogenetic --i-phylogeny %s ' % tree
    elif metric not in ALPHA_PHYLOGENETIC_METRICS and tree is None:
        cmd = 'qiime diversity alpha '
    else:
        return False, None, ('Phylogenetic metric %s selected but no tree '
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

    ainfo = [ArtifactInfo('Alpha vectors', 'alpha_vector',
                          [(ffp, 'plain_text')])]
    return True, ainfo, ""


def alpha_correlation(qclient, job_id, parameters, out_dir):
    """generate alpha correlation calculations

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values for alpha correlation
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
    artifact_id = parameters['Alpha vectors']
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    dm_fp = artifact_info['files']['plain_text'][0]
    dm_qza = join(out_dir, 'q2-alpha-diversity.qza')
    analysis_id = artifact_info['analysis']
    metadata = qclient.get(
        "/qiita_db/analysis/%s/metadata/" % str(analysis_id))
    metadata = pd.DataFrame.from_dict(metadata, orient='index')
    metadata_fp = join(out_dir, 'metadata.txt')
    metadata.to_csv(metadata_fp, sep='\t')
    p_method = ALPHA_CORRELATION_METHODS[parameters['Correlation method']]
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

    ainfo = [ArtifactInfo('Alpha correlation visualization',
                          'q2_visualization', [(o_visualization, 'qzv')])]
    return True, ainfo, ""


def taxa_barplot(qclient, job_id, parameters, out_dir):
    """Generate taxa barplot calculations

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values for taxa barplot
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    boolean, list, str
        The results of the job
    """
    out_dir = join(out_dir, 'taxa_barplot')
    if not exists(out_dir):
        mkdir(out_dir)

    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    artifact_id = int(parameters['BIOM table'])
    artifact_info = qclient.get("/qiita_db/artifacts/%d/" % artifact_id)
    analysis_id = artifact_info['analysis']
    metadata = qclient.get(
        "/qiita_db/analysis/%s/metadata/" % str(analysis_id))
    metadata = pd.DataFrame.from_dict(metadata, orient='index')
    metadata_fp = join(out_dir, 'metadata.txt')
    metadata.to_csv(metadata_fp, sep='\t')

    biom_qza = join(out_dir, 'q2-biom.qza')
    taxonomy_txt = join(out_dir, 'q2-taxonomy.txt')
    taxonomy_qza = join(out_dir, 'q2-taxonomy.qza')
    taxa_plot_qzv = join(out_dir, 'taxa-barplot.qzv')

    # getting the biom table so we can check for taxonomies
    biom_fp = artifact_info['files']['biom'][0]
    bt = load_table(biom_fp)
    with open(taxonomy_txt, 'w') as fp:
        fp.write('Feature ID\tTaxon\n')
        for otu_id in bt.ids('observation'):
            tax = bt.metadata(id=otu_id, axis='observation')
            if tax is None:
                error_msg = ("biom table doesn't have taxonomy")
                return False, None, error_msg
            taxonomy = '; '.join(tax['taxonomy'])
            fp.write("%s\t%s\n" % (otu_id, taxonomy))

    qclient.update_job_step(
        job_id, "Step 2 of 4: Converting Qiita artifacts to Q2 artifact: BIOM")
    cmd = ('qiime tools import --input-path %s --output-path %s '
           '--type "FeatureTable[Frequency]' % (biom_fp, biom_qza))

    counts = list(map(sum, bt.iter_data()))
    if min(counts) == max(counts):
        cmd += " % Properties(['uniform-sampling'])\""
    else:
        cmd += '"'
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error converting biom:\nStd out: %s\nStd err: %s"
                     % (std_out, std_err))
        return False, None, error_msg

    qclient.update_job_step(job_id, "Step 3 of 4: Converting Qiita artifacts "
                                    "to Q2 artifact: Taxonomy")
    cmd = ('qiime tools import --input-path %s --output-path %s '
           '--type "FeatureData[Taxonomy]"' % (taxonomy_txt, taxonomy_qza))
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error converting taxonomy:\nStd out: %s\n"
                     "Std err: %s" % (std_out, std_err))
        return False, None, error_msg

    qclient.update_job_step(job_id, "Step 4 of 4: Generating summary")
    cmd = ('qiime taxa barplot --i-table %s --i-taxonomy %s '
           '--m-metadata-file %s --o-visualization %s' % (
               biom_qza, taxonomy_qza, metadata_fp, taxa_plot_qzv))
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error generating taxonomy summary:\nStd out: %s\n"
                     "Std err: %s" % (std_out, std_err))
        return False, None, error_msg

    ainfo = [ArtifactInfo('Taxa summaries visualization', 'q2_visualization',
                          [(taxa_plot_qzv, 'qzv')])]
    return True, ainfo, ""


def filter_samples(qclient, job_id, parameters, out_dir):
    """Filter samples from a table

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values for filter samples
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    boolean, list, str
        The results of the job
    """
    out_dir = join(out_dir, 'filter_samples')
    if not exists(out_dir):
        mkdir(out_dir)

    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    artifact_id = int(parameters['BIOM table'])
    p_max_frequency = int(
        parameters['Maximum feature frequency across samples'])
    p_max_features = int(parameters['Maximum features per sample'])
    p_min_frequency = int(
        parameters['Minimum feature frequency across samples'])
    p_min_features = int(parameters['Minimum features per sample'])
    p_where = parameters['SQLite WHERE-clause']

    artifact_info = qclient.get("/qiita_db/artifacts/%d/" % artifact_id)
    analysis_id = artifact_info['analysis']
    metadata = qclient.get(
        "/qiita_db/analysis/%s/metadata/" % str(analysis_id))
    metadata = pd.DataFrame.from_dict(metadata, orient='index')
    metadata_fp = join(out_dir, 'metadata.txt')
    metadata.to_csv(metadata_fp, sep='\t')

    # getting just the biom file, [0] it should be only one
    biom_ifp = artifact_info['files']['biom'][0]
    biom_ofp = join(out_dir, 'biom.qza')

    qclient.update_job_step(
        job_id, "Step 2 of 4: Converting Qiita artifacts to Q2 artifact")
    # converting biom
    cmd = ('qiime tools import --input-path %s --output-path %s '
           '--type "FeatureTable[Frequency]' % (biom_ifp, biom_ofp))
    b = load_table(biom_ifp)
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

    qclient.update_job_step(job_id, "Step 3 of 4: Filtering")
    filter_ofp = join(out_dir, 'biom_filtered.qza')
    if 'Exclude ids selected by where parameter':
        exclude_ids = '--p-no-exclude-ids'
    else:
        exclude_ids = '--p-exclude-ids'
    cmd = ('qiime feature-table filter-samples --m-metadata-file %s '
           '--o-filtered-table %s --p-max-frequency %d --p-max-features %d '
           '--p-min-frequency %d --p-min-features %d --i-table %s %s' % (
               metadata_fp, filter_ofp, p_max_frequency, p_max_features,
               p_min_frequency, p_min_features, biom_ofp, exclude_ids))
    if p_where != '':
        cmd += ' --p-where "%s"' % p_where
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error in filtering samples in biom\nStd out: %s\n"
                     "Std err: %s" % (std_out, std_err))
        return False, None, error_msg

    qclient.update_job_step(
        job_id, "Step 4 of 4: Converting Q2 to Qiita artifacts")
    fdir = join(out_dir, 'filter_samples')
    ffp = join(fdir, 'feature-table.biom')
    cmd = "qiime tools export --output-dir %s %s" % (fdir, filter_ofp)
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error in Q2 -> Qiita conversion:\nStd out: "
                     "%s\nStd err: %s" % (std_out, std_err))
        return False, None, error_msg

    # After calling Qiime2, the taxonomy has been dropped from the BIOM table
    # Re-add here
    orig = load_table(biom_ifp)
    res = load_table(ffp)

    metadata = {i: orig.metadata(i, axis='observation')
                for i in res.ids(axis='observation')}
    res.add_metadata(metadata, axis='observation')

    res_fp = join(out_dir, 'filtered.biom')
    with biom_open(res_fp, 'w') as bf:
        res.to_hdf5(bf, "Qiita's Qiime2 plugin")

    ainfo = [ArtifactInfo('Filtered table', 'BIOM', [(res_fp, 'biom')])]
    return True, ainfo, ""


def emperor(qclient, job_id, parameters, out_dir):
    """generate emperor plot calculations

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values for pcoa
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    boolean, list, str
        The results of the job
    """
    out_dir = join(out_dir, 'emperor')
    if not exists(out_dir):
        mkdir(out_dir)

    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    artifact_id = parameters['Ordination results']
    p_custom_axis = parameters['Custom axis']
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    pcoa_fp = artifact_info['files']['plain_text'][0]

    analysis_id = artifact_info['analysis']
    metadata = qclient.get(
        "/qiita_db/analysis/%s/metadata/" % str(analysis_id))
    metadata = pd.DataFrame.from_dict(metadata, orient='index')
    metadata_fp = join(out_dir, 'metadata.txt')
    metadata.to_csv(metadata_fp, sep='\t')

    pcoa_qza = join(out_dir, 'q2-pcoa.qza')
    emperor_qzv = join(out_dir, 'q2-emperor.qzv')

    qclient.update_job_step(
        job_id, "Step 2 of 4: Converting Qiita artifacts to Q2 artifact")
    cmd = ('qiime tools import --input-path %s --output-path %s '
           '--type "PCoAResults"' % (pcoa_fp, pcoa_qza))
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error converting distance matrix:\nStd out: %s\n"
                     "Std err: %s" % (std_out, std_err))
        return False, None, error_msg

    qclient.update_job_step(
        job_id, "Step 3 of 4: Generating Emperor plot")

    cmd = ('qiime emperor plot --i-pcoa %s --o-visualization %s '
           '--m-metadata-file %s' % (pcoa_qza, emperor_qzv, metadata_fp))
    if p_custom_axis is not None and p_custom_axis not in ['None', '']:
        cmd += ' --p-custom-axis "%s"' % p_custom_axis

    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error in PCoA\nStd out: %s\nStd err: %s"
                     % (std_out, std_err))
        return False, None, error_msg

    ainfo = [ArtifactInfo('Emperor visualization', 'q2_visualization',
                          [(emperor_qzv, 'qzv')])]
    return True, ainfo, ""


def beta_group_significance(qclient, job_id, parameters, out_dir):
    """generate beta group significance calculations

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values for beta correlation
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    boolean, list, str
        The results of the job
    """
    out_dir = join(out_dir, 'beta_group_significance')
    if not exists(out_dir):
        mkdir(out_dir)

    qclient.update_job_step(job_id, "Step 1 of 3: Collecting information")
    artifact_id = parameters['Distance matrix']
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    dm_fp = artifact_info['files']['plain_text'][0]
    dm_qza = join(out_dir, 'q2-distance.qza')
    analysis_id = artifact_info['analysis']
    metadata = qclient.get(
        "/qiita_db/analysis/%s/metadata/" % str(analysis_id))
    metadata = pd.DataFrame.from_dict(metadata, orient='index')
    metadata_fp = join(out_dir, 'metadata.txt')
    metadata.to_csv(metadata_fp, sep='\t')
    m_metadata_category = parameters['Metadata category']
    p_method = BETA_GROUP_SIG_METHODS[parameters['Method']]
    p_permutations = parameters['Number of permutations']
    p_pairwise = BETA_GROUP_SIG_TYPE[parameters['Comparison type']]
    o_visualization = join(out_dir, 'beta_group_significance.qzv')

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
        job_id, "Step 3 of 3: Calculating beta group significance")
    cmd = ('qiime diversity beta-group-significance --i-distance-matrix %s '
           '--m-metadata-file %s --m-metadata-category %s --p-method %s '
           '--p-permutations %s --o-visualization %s --%s' % (
               dm_qza, metadata_fp, m_metadata_category, p_method,
               p_permutations, o_visualization, p_pairwise))
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error in beta group significance\nStd out: %s\n"
                     "Std err: %s" % (std_out, std_err))
        return False, None, error_msg

    ainfo = [ArtifactInfo('Beta group significance visualization',
                          'q2_visualization',
                          [(o_visualization, 'qzv')])]
    return True, ainfo, ""
