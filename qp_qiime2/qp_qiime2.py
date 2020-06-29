# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os import mkdir, listdir, chmod
from os.path import join, exists, basename
from shutil import copyfile

from biom import load_table
from biom.util import biom_open

from qiita_client import ArtifactInfo

import qiime2
import pandas as pd


Q2_ALLOWED_PLUGINS = [
    'taxa', 'sample-classifier', 'composition', 'phylogeny', 'feature-table',
    'gneiss', 'diversity', 'longitudinal', 'emperor'
]

QIITA_Q2_SEMANTIC_TYPE = {
    'BIOM': {
        'name': 'FeatureTable',
        'expression': ['Frequency', 'RelativeFrequency', 'PresenceAbsence']},
    'distance_matrix': {
        'name': 'DistanceMatrix',
        'expression': []},
    'ordination_results': {
        'name': 'PCoAResults',
        'expression': []},
    'q2_visualization': {
        'name': 'Visualization',
        'expression': []},
    'alpha_vector': {
        'name': 'SampleData',
        'expression': ['AlphaDiversity']},
    'phylogeny': {
        'name': 'Phylogeny',
        'expression': ['Rooted']},
    'FeatureData[Taxonomy]':  {
        'name': 'FeatureData',
        'expression': ['Taxonomy']},
}

# for simplicity we are going to invert QIITA_Q2_SEMANTIC_TYPE so we can
# search by key or value without having to do this operation several times
Q2_QIITA_SEMANTIC_TYPE = {
    y['name']: x for x, y in QIITA_Q2_SEMANTIC_TYPE.items()}

PRIMITIVE_TYPES = {
    'Int': 'integer',
    'Str': 'string',
    'Choices': 'choice',
    'Set': 'choice',
    'List': 'choice',
    'Metadata': 'string',
    'MetadataColumn': 'string',
    'Float': 'float',
    'Bool': 'boolean',
}

ALPHA_DIVERSITY_METRICS_PHYLOGENETIC = {
    "Faith's Phylogenetic Diversity": "faith_pd",
}

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

BETA_DIVERSITY_METRICS = {
    "Aitchison distance": "aitchison",
    "Bray-Curtis dissimilarity": "braycurtis",
    "Canberra distance": "canberra",
    "Canberra distance in Adkins form": "canberra_adkins",
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
    "Jensen-Shannon metric": "jensenshannon",
    "Yule index": "yule"}

BETA_DIVERSITY_METRICS_PHYLOGENETIC = {
    "Unweighted UniFrac": "unweighted_unifrac",
    'Weighted UniFrac': 'weighted_unifrac',
    "Weighted normalized UniFrac": "weighted_normalized_unifrac",
    "Generalized UniFrac": "generalized_unifrac"
}

CORRELATION_METHODS = {
    "Spearman": "spearman",
    "Pearson": "pearson"}

BETA_GROUP_SIG_METHODS = {
    "PERMDISP": "permdisp",
    "PERMANOVA": "permanova",
    "ANOSIM": "anosim"}

RENAME_COMMANDS = {
    ('alpha', 'metric'): ALPHA_DIVERSITY_METRICS,
    ('beta', 'metric'): BETA_DIVERSITY_METRICS,
    ('alpha_phylogenetic', 'metric'): ALPHA_DIVERSITY_METRICS_PHYLOGENETIC,
    ('beta_phylogenetic', 'metric'): BETA_DIVERSITY_METRICS_PHYLOGENETIC,
    ('alpha_rarefaction', 'metrics'): {
        **ALPHA_DIVERSITY_METRICS, **ALPHA_DIVERSITY_METRICS_PHYLOGENETIC},
    ('beta_rarefaction', 'metric'): {
        **BETA_DIVERSITY_METRICS, **BETA_DIVERSITY_METRICS_PHYLOGENETIC},
    ('beta_rarefaction', 'correlation_method'): CORRELATION_METHODS,
    ('beta_correlation', 'method'): CORRELATION_METHODS,
    ('alpha_correlation', 'method'): CORRELATION_METHODS,
    ('beta_group_significance', 'method'): BETA_GROUP_SIG_METHODS,
}


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
    # going to separate in q2params and q2inputs as the inputs are going to
    # need to be retrieved from qiita and converted to qza
    label = 'qp-hide-param'
    label_len = len(label)
    q2params = {}
    q2inputs = {}
    method_inputs = method.signature.inputs.copy()
    method_params = method.signature.parameters.copy()
    artifact_id = None
    analysis_id = None
    biom_fp = None
    tree_fp = None
    tree_fp_check = False
    for k in list(parameters):
        if k in parameters and k.startswith(label):
            key = parameters.pop(k)
            val = parameters.pop(k[label_len:])
            if key in method_inputs.keys():
                if key == 'phylogeny':
                    if val == '':
                        continue
                    # there is a chance that we parse/loop over the phylogeny
                    # option before the artifact so tree_fp will still be
                    # None; thus we will need to check this after we are done
                    # with this loop
                    if val == 'Artifact tree, if exists':
                        tree_fp_check = True
                    fpath = val
                    qiita_name = QIITA_Q2_SEMANTIC_TYPE[key]
                    if qiita_name['expression']:
                        # for these cases we need an expresion so for
                        # simplicity using the first one [0]
                        artifact_method = '%s[%s]' % (
                            qiita_name['name'], qiita_name['expression'][0])
                elif key == 'classifier':
                    fpath = val
                    artifact_method = None
                    k = key
                else:
                    # this is going to be an artifact so let's collect the
                    # filepath here, this will also allow us to collect the
                    # analysis_id
                    artifact_id = val
                    ainfo = qclient.get(
                        "/qiita_db/artifacts/%s/" % artifact_id)
                    if ainfo['analysis'] is None:
                        msg = ('Artifact "%s" is not an analysis '
                               'artifact.' % val)
                        return False, None, msg
                    analysis_id = ainfo['analysis']
                    dt = method_inputs[key].qiime_type.to_ast()['name']
                    if 'qza' not in ainfo['files']:
                        # at this stage in qiita we only have 2 types of
                        # artifacts: biom / plain_text
                        if Q2_QIITA_SEMANTIC_TYPE[dt] == 'BIOM':
                            fpath = ainfo['files']['biom'][0]
                            biom_fp = fpath
                        else:
                            fpath = ainfo['files']['plain_text'][0]
                    else:
                        fpath = ainfo['files']['qza'][0]
                    # if it's a BIOM and there is a plain_text is the
                    # result of the archive at this stage: a tree
                    if Q2_QIITA_SEMANTIC_TYPE[dt] == 'BIOM':
                        if 'plain_text' in ainfo['files']:
                            tree_fp = ainfo['files']['plain_text'][0]
                    if biom_fp is None and 'biom' in ainfo['files']:
                        biom_fp = ainfo['files']['biom'][0]

                    q2artifact_name = Q2_QIITA_SEMANTIC_TYPE[
                        method_inputs[key].qiime_type.to_ast()['name']]
                    qiita_name = QIITA_Q2_SEMANTIC_TYPE[q2artifact_name]
                    if qiita_name['expression']:
                        # for these cases we need an expresion so for
                        # simplicity using the first one [0]
                        artifact_method = '%s[%s]' % (
                            qiita_name['name'], qiita_name['expression'][0])
                    else:
                        artifact_method = qiita_name['name']

                q2inputs[key] = (fpath, artifact_method)
            elif key == 'qp-hide-metadata-field':
                if val == '':
                    msg = ("Error: You didn't write a metadata field in "
                           "'%s'" % k[label_len:])
                    return False, None, msg
                q2inputs['metadata'] = (val, val)
            else:
                if val in ('', 'None'):
                    continue

                mkey = method_params[key]
                # users can only select one value so if the view_type is set
                # we will not convert
                if mkey.view_type not in (set, str):
                    if mkey.view_type is not mkey.NOVALUE:
                        val = mkey.view_type(val)
                    elif val not in mkey.qiime_type:
                        val = mkey.qiime_type.decode(val)

                # let's bring back the original name of these parameters
                value_pair = (q2method, key)
                if (q2plugin == 'diversity' and value_pair in RENAME_COMMANDS):
                    val = RENAME_COMMANDS[value_pair][val]

                # if the view_type is set we need to convert to an actual
                # set
                if mkey.view_type is set:
                    val = {val}
                q2params[key] = val
        elif k in ('qp-hide-metadata', 'qp-hide-FeatureData[Taxonomy]'):
            # remember, if we need metadata, we will always have
            # qp-hide-metadata and optionaly we will have
            # qp-hide-metadata-field
            key = parameters.pop(k)
            q2inputs[key] = ('', '')

    # if 'metadata' is in q2inputs but 'where' exist and is empty in q2params,
    # remove the parameter metadata
    # NOTE: AFAIK there is no way to differentiate between sample and prep
    #       metadata in Q2 so the need to remove for filter_features

    if ('filter_features' == q2method or (
            'metadata' in q2inputs and ('where' in q2params and
                                        not q2params['where']))):
        q2inputs.pop('metadata')

    # if we are here, we need to use the internal tree from the artifact
    if tree_fp_check:
        q2inputs['phylogeny'] = (tree_fp, q2inputs['phylogeny'][1])

    # let's process/import inputs
    qclient.update_job_step(
        job_id, "Step 2 of 4: Converting Qiita artifacts to Q2 artifact")
    for k, (fpath, dt) in q2inputs.items():
        if k in ('metadata', 'sample_metadata'):
            metadata = qclient.get(
                "/qiita_db/analysis/%s/metadata/" % str(analysis_id))
            metadata = pd.DataFrame.from_dict(metadata, orient='index')
            # the reason we need to save and load the mapping file is
            # so Qiime2 assings the expected data types to the columns
            metadata_fp = join(out_dir, 'metadata.txt')
            metadata.to_csv(metadata_fp, index_label='#SampleID', na_rep='',
                            sep='\t', encoding='utf-8')
            q2Metadata = qiime2.Metadata.load(metadata_fp)
            if fpath:
                q2params[k] = q2Metadata.get_column(fpath)
            else:
                q2params[k] = q2Metadata
        elif k == 'FeatureData[Taxonomy]':
            try:
                qza = qiime2.Artifact.import_data(
                    'FeatureData[Taxonomy]', biom_fp, 'BIOMV210Format')
            except Exception:
                return False, None, ('Error generating taxonomy. Are you '
                                     'sure this artifact has taxonomy?')
            q2params['taxonomy'] = qza
        else:
            if not fpath.endswith('.qza'):
                try:
                    qza = qiime2.Artifact.import_data(dt, fpath)
                except Exception as e:
                    return False, None, 'Error converting "%s": %s' % (
                        str(dt), str(e))
            else:
                qza = qiime2.Artifact.load(fpath)
            q2params[k] = qza

    # if feature_classifier and classify_sklearn we need to transform the
    # input data to sequences
    if q2plugin == 'feature-classifier' and q2method == 'classify_sklearn':
        ainfo = qclient.get("/qiita_db/artifacts/%s/" %
                            parameters['The feature data to be classified.'])
        biom_fp = ainfo['files']['biom'][0]
        plain_text_fp = None
        if 'plain_text' in ainfo['files']:
            plain_text_fp = ainfo['files']['plain_text'][0]
        biom_table = load_table(biom_fp)
        fna_fp = join(out_dir, 'sequences.fna')
        with open(fna_fp, 'w') as f:
            for _id in biom_table.ids(axis='observation'):
                f.write('>{0}\n{0}\n'.format(_id))
        try:
            q2params['reads'] = qiime2.Artifact.import_data(
                'FeatureData[Sequence]', fna_fp)
        except (ValueError, qiime2.core.exceptions.ValidationError) as e:
            msg = str(e)
            if '(does not match IUPAC characters for a DNA sequence)' in msg:
                msg = ('Table IDs are not sequences, please confirm that this '
                       'is not a closed reference table?')
            return False, None, 'Error converting "%s": %s' % (
                'Input Table', msg)

    qclient.update_job_step(
        job_id, "Step 3 of 4: Running '%s %s'" % (q2plugin, q2method))
    try:
        results = method(**q2params)
    except Exception as e:
        return False, None, 'Error running: %s' % str(e)

    qclient.update_job_step(job_id, "Step 4 of 4: Processing results")
    out_info = []

    # if feature_classifier and classify_sklearn we need to add the taxonomy
    # to the original table and generate the new artifact
    if q2plugin == 'feature-classifier' and q2method == 'classify_sklearn':
        new_biom = join(out_dir, 'feature-table-with-taxonomy.biom')
        new_qza = join(out_dir, 'feature-table-with-taxonomy.qza')
        df = results[0].view(pd.DataFrame)
        df.rename(columns={'Taxon': 'taxonomy'}, inplace=True)
        df['taxonomy'] = [[y.strip() for y in x]
                          for x in df['taxonomy'].str.split(';')]
        biom_table.add_metadata(df.to_dict(orient='index'), axis='observation')
        with biom_open(new_biom, 'w') as bf:
            biom_table.to_hdf5(bf, 'Generated in Qiita')

        qza = qiime2.Artifact.import_data(
            'FeatureTable[Frequency]', new_biom, 'BIOMV210Format')
        qza.save(new_qza)
        ftc_fps = [(new_biom, 'biom'), (new_qza, 'qza')]
        if plain_text_fp is not None:
            # if we enter here, it means that the input artifact had a tree
            # (saved as plain_text); thus, we need to make sure we make a copy
            # so we don't move the original file
            bn = basename(plain_text_fp)
            new_tree_fp = join(out_dir, bn)
            copyfile(ainfo['files']['plain_text'][0], new_tree_fp)
            ftc_fps.append((new_tree_fp, 'plain_text'))
        out_info.append(ArtifactInfo(
            'Feature Table with Classification', 'BIOM', ftc_fps))

    for aname, q2artifact in zip(results._fields, results):
        aout = join(out_dir, aname)
        if isinstance(q2artifact, qiime2.Visualization):
            qzv_fp = q2artifact.save(aout)
            out_info.append(
                ArtifactInfo(aname, 'q2_visualization', [(qzv_fp, 'qzv')]))
        else:
            qza_fp = q2artifact.save(aout + '.qza')
            q2artifact.export_data(output_dir=aout)
            files = listdir(aout)
            if len(files) != 1:
                msg = ('Error processing results: There are some unexpected '
                       'files: "%s"' % ', '.join(files))
                return False, None, msg
            fp = join(aout, files[0])
            # making sure the newly created file comes with the correct
            # permissions for nginx
            chmod(fp, 0o664)

            if (q2artifact.type.name == 'FeatureTable'):
                # Re-add the observation metadata if exists in the input and if
                # not one of the plugin/methods that actually changes that
                # information
                if biom_fp is not None and (q2plugin, q2method) not in [
                        ('taxa', 'collapse')]:
                    fin = load_table(biom_fp)
                    fout = load_table(fp)

                    # making sure that the resulting biom is not empty
                    if fout.shape == (0, 0):
                        msg = ('The resulting table is empty, please review '
                               'your parameters')
                        return False, None, msg

                    metadata = {
                        i: fin.metadata(i, axis='observation')
                        for i in fout.ids(axis='observation')}
                    fout.add_metadata(metadata, axis='observation')
                    with biom_open(fp, 'w') as bf:
                        fout.to_hdf5(bf, "Qiita's Qiime2 plugin with "
                                     "observation metadata")

                # if there is a tree, let's copy it and then add it to
                # the new artifact
                if tree_fp is not None:
                    bn = basename(tree_fp)
                    new_tree_fp = join(
                        out_dir, aout, 'from_%s_%s' % (artifact_id, bn))
                    copyfile(tree_fp, new_tree_fp)
                    ai = ArtifactInfo(aname, 'BIOM', [
                        (fp, 'biom'),
                        (new_tree_fp, 'plain_text'),
                        (qza_fp, 'qza')])
                else:
                    ai = ArtifactInfo(
                        aname, 'BIOM', [(fp, 'biom'), (qza_fp, 'qza')])

            else:
                atype = Q2_QIITA_SEMANTIC_TYPE[q2artifact.type.name]
                ai = ArtifactInfo(
                    aname, atype, [(fp, 'plain_text'), (qza_fp, 'qza')])
            out_info.append(ai)

    return True, out_info, ""
