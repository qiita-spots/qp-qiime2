# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os import mkdir, listdir
from os.path import join, exists

from qiita_client import ArtifactInfo

import qiime2
import pandas as pd


Q2_ALLOWED_PLUGINS = [
    'taxa', 'sample-classifier', 'composition', 'phylogeny', 'feature-table',
    'gneiss', 'diversity', 'longitudinal', 'emperor'
]

QIITA_Q2_SEMANTIC_TYPE = {
    'BIOM-F': qiime2.sdk.util.parse_type('FeatureTable[Frequency]'),
    'BIOM-RF': qiime2.sdk.util.parse_type('FeatureTable[RelativeFrequency]'),
    'BIOM-PA': qiime2.sdk.util.parse_type('FeatureTable[PresenceAbsence]'),
    'distance_matrix': qiime2.sdk.util.parse_type('DistanceMatrix'),
    'ordination_results': qiime2.sdk.util.parse_type('PCoAResults'),
    'q2_visualization': qiime2.sdk.util.parse_type('Visualization'),
    'alpha_vector': qiime2.sdk.util.parse_type('SampleData[AlphaDiversity]'),
    'phylogenetic_distance_matrix': qiime2.sdk.util.parse_type(
        "DistanceMatrix % Properties(['phylogenetic'])"),
    'phylogenetic_alpha_vector': qiime2.sdk.util.parse_type(
        "SampleData[AlphaDiversity] % Properties(['phylogenetic'])"),
    'phylogeny': qiime2.sdk.util.parse_type('Phylogeny[Rooted]'),
    'taxonomy': qiime2.sdk.util.parse_type('FeatureData[Taxonomy]')
}

# for simplicity we are going to invert QIITA_Q2_SEMANTIC_TYPE so we can
# search by key or value without having to do this operation several times
Q2_QIITA_SEMANTIC_TYPE = {y: x for x, y in QIITA_Q2_SEMANTIC_TYPE.items()}

PRIMITIVE_TYPES = {
    qiime2.core.type.primitive._Int: 'integer',
    qiime2.core.type.primitive._Bool: 'boolean',
    qiime2.core.type.primitive._Float: 'float',
    qiime2.core.type.primitive._Str: 'string',
    qiime2.core.type.collection._CollectionPrimitive: 'choice',
    # this is a mapping file
    qiime2.core.type.primitive._Metadata: 'mapping',
    # this is a column in the mapping file
    qiime2.core.type.primitive._MetadataColumnExpression: 'string',
}

STATE_UNIFRAC_METRICS = {
    "unweighted UniFrac": "unweighted_unifrac",
    "weighted normalized UniFrac": "weighted_normalized_unifrac",
    "weighted unnormalized UniFrac": "weighted_unifrac",
    "generalized UniFrac": "generalized_unifrac",
}

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
    "Chebyshev distance": "chebyshev",
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


def call_qiime2(qclient, job_id, parameters, out_dir):
    """helper method to call Qiime2

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to process
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    boolean, list, str
        The results of the job
    """
    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    q2plugin = parameters.pop('qp-hide-plugin')
    q2method = parameters.pop('qp-hide-method').replace('-', '_')
    pm = qiime2.sdk.PluginManager()
    method = pm.plugins[q2plugin].actions[q2method]

    out_dir = join(out_dir, q2method)

    # making sure that we always start with an empty folder
    if not exists(out_dir):
        mkdir(out_dir)

    # let's generate the parameters, first remove the hidden parameters. We are
    # gonna separate in q2params and q2inputs as the inputs are going to need
    # to be retrieved from qiita and converted to qza
    label = 'qp-hide-param'
    label_len = len(label)
    q2params = {}
    q2inputs = {}
    method_inputs = method.signature.inputs.copy()
    method_params = method.signature.parameters.copy()
    analysis_id = None
    for k in list(parameters):
        if k in parameters and k.startswith(label):
            key = parameters.pop(k)
            val = parameters.pop(k[label_len:])
            if key in method_inputs.keys():
                if key in ('phylogeny'):
                    # ******** ToDo **********
                    # here we need to add taxonomy
                    # ************************
                    fpath = val
                    artifact_method = QIITA_Q2_SEMANTIC_TYPE[key]
                else:
                    # this is going to be an artifact so let's collect the
                    # filepath here, this will also allow us to collect the
                    # analysis_id
                    ainfo = qclient.get("/qiita_db/artifacts/%s/" % val)
                    if ainfo['analysis'] is None:
                        msg = ('Artifact "%s" is not an analysis '
                               'artifact.' % val)
                        return False, None, msg
                    analysis_id = ainfo['analysis']
                    # at this stage in qiita we only have 2 types of artifacts:
                    # biom / plain_text
                    dt = method_inputs[key].qiime_type
                    fp = 'plain_text'
                    if Q2_QIITA_SEMANTIC_TYPE[dt].startswith('BIOM'):
                        fp = 'biom'
                    fpath = ainfo['files'][fp][0]
                    artifact_method = method_inputs[key].qiime_type
                q2inputs[key] = (fpath, artifact_method)
            else:
                if val in ('', 'None'):
                    continue
                val = method_params[key].view_type(val)
                q2params[key] = val
        elif k == 'qp-hide-metadata':
            # remember, if we need metadata, we will always have
            # qp-hide-metadata and optionaly we will have
            # qp-hide-metadata-field
            key = parameters.pop(k)
            if 'qp-hide-metadata-field' in parameters:
                val = parameters.pop(k[label_len:])
                if val == '':
                    continue
                q2params[key] = val
            q2inputs[key] = ('', '')

    # let's process/import inputs
    qclient.update_job_step(
        job_id, "Step 2 of 4: Converting Qiita artifacts to Q2 artifact")
    for k, (fpath, dt) in q2inputs.items():
        if k == 'metadata':
            metadata = qclient.get(
                "/qiita_db/analysis/%s/metadata/" % str(analysis_id))
            metadata = pd.DataFrame.from_dict(metadata, orient='index')
            # the reason we need to save and load the mapping file is
            # so Qiime2 assings the expected data types to the columns
            metadata_fp = join(out_dir, 'metadata.txt')
            metadata.to_csv(metadata_fp, index_label='#SampleID', na_rep='',
                            sep='\t', encoding='utf-8')
            q2params[k] = qiime2.Metadata.load(metadata_fp)
        else:
            try:
                qza = qiime2.Artifact.import_data(dt, fpath)
            except Exception as e:
                return False, None, 'Error converting "%s": %s' % (
                    str(dt), str(e))
            q2params[k] = qza

    qclient.update_job_step(
        job_id, "Step 3 of 4: Running '%s %s'" % (q2plugin, q2method))
    try:
        results = method(**q2params)
    except Exception as e:
        return False, None, 'Error running: %s' % str(e)

    qclient.update_job_step(job_id, "Step 4 of 4: Processing results")
    ainfo = []
    for aname, q2artifact in zip(results._fields, results):
        aout = join(out_dir, aname)
        q2artifact.export_data(output_dir=aout)
        if isinstance(q2artifact, qiime2.Visualization):
            ainfo.append(
                ArtifactInfo(aname, 'q2_visualization', [(aout, 'qzv')]))
        else:
            files = listdir(aout)
            if len(files) != 1:
                msg = ('Error processing results: There are some unexpected '
                       'files: "%s"' % ', '.join(files))
                return False, None, msg
            fp = join(aout, files[0])

            if q2artifact.type.name == 'FeatureTable':
                atype = 'BIOM'
                ftype = 'biom'
            else:
                atype = Q2_QIITA_SEMANTIC_TYPE[qza.type]
                ftype = 'plain_text'
            ainfo.append(ArtifactInfo(aname, atype, [(fp, ftype)]))

    return True, ainfo, ""