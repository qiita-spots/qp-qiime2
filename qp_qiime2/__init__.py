# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from json import dumps
from os import environ
from os.path import join
from glob import glob

from qiita_client import QiitaPlugin, QiitaCommand

from .qp_qiime2 import (
    QIITA_Q2_SEMANTIC_TYPE, Q2_QIITA_SEMANTIC_TYPE, Q2_ALLOWED_PLUGINS,
    PRIMITIVE_TYPES, call_qiime2, RENAME_COMMANDS, NOT_VALID_OUTPUTS,
    Q2_EXTRA_COMMANDS)
from .util import get_qiime2_type_name_and_predicate
from qiime2 import __version__ as qiime2_version
from qiime2.sdk.util import actions_by_input_type
from qiime2.sdk import PluginManager


# Initialize the qiita_plugin
plugin = QiitaPlugin('qiime2', qiime2_version, 'QIIME 2')

# The extra commands require a folder where all the pre-calculated databases
# exist, which is set up via a ENV variable, if not present we should raise
# an error
qp_qiime2_dbs = environ.get('QP_QIIME2_DBS')
if qp_qiime2_dbs is None:
    raise ValueError("Missing ENV var QP_QIIME2_DBS, please set.")
qp_qiime2_dbs = glob(join(qp_qiime2_dbs, '*.qza'))
if len(qp_qiime2_dbs) < 1:
    raise ValueError(
        "ENV QP_QIIME2_DBS points to a folder without QZA files, please set.")

# Similar to the precomputed databases, we also need a folder containing all
# the valid per-sequence filtering artifacts
qp_filtering_qza = environ.get('QP_QIIME2_FILTER_QZA')
if qp_filtering_qza is None:
    raise ValueError("Missing ENV var QP_QIIME2_FILTER_QZA, please set.")
qp_filtering_qza = glob(join(qp_filtering_qza, '*.qza'))
if len(qp_filtering_qza) < 1:
    raise ValueError("ENV QP_QIIME2_FILTER_QZA points to a folder without "
                     "QZA files, please set.")

# PLEASE READ:
# There are 2 main steps:
# 1. We are going to loop on the Q2_EXTRA_COMMANDS to check which extra
# commands we want to add, extra from those given in the next step
# 2. We are going loop over QIITA_Q2_SEMANTIC_TYPE (lookup table) so we can
# retrieve the q2plugin and their methods that work with that given Q2/Qiita
# semantic type. Then we will ignore any plugin not in Q2_ALLOWED_PLUGINS so
# we avoid adding plugins that we don't want; like deblur or dada2. Finally,
# we are going to loop over the different inputs, outputs and parameters from
# Q2 and convert them to QIITA's req_params, opt_params and outputs.
#
# Note that Qiita users like to have descriptions of the paramters
# (q2-description) vs. the parameter itself (q2-parameter) so to allow this
# we are going to store each parameter twice: one in the
# opt_params[q2-description]: value; and
# req_params['qp-hide-param' + q2-description]: q2-parameter

pm = PluginManager()
methods_to_add = []
for plugin_name, method_name in Q2_EXTRA_COMMANDS:
    q2plugin = pm.plugins[plugin_name]
    m = q2plugin.methods[method_name]
    methods_to_add.append((q2plugin, m))

for qiita_artifact, q2_artifacts in QIITA_Q2_SEMANTIC_TYPE.items():
    if q2_artifacts['expression']:
        actions = [a for e in q2_artifacts['expression']
                   for a in actions_by_input_type('%s[%s]' % (
                       q2_artifacts['name'], e))]
    else:
        actions = actions_by_input_type(q2_artifacts['name'])

    for q2plugin, methods in actions:
        # note that the qiita_artifact are strings not objects
        if qiita_artifact.startswith('BIOM'):
            qiita_artifact = 'BIOM'

        if q2plugin.name not in Q2_ALLOWED_PLUGINS:
            # As of qiime2-2022.11 this filters out:
            # alignment
            # deblur
            # diversity-lib
            # feature-classifier
            # fragment-insertion
            # quality-control
            # sourcetracker2
            # vsearch
            continue

        for m in methods:
            # after review of qiime2-2019.4 we decided to not add these methods
            if (q2plugin.name, m.id) not in [('feature-table', 'group'),
                                             ('feature-table', 'filter_seqs')]:
                methods_to_add.append((q2plugin, m))

# make sure we have seen all expected plugins
q2_expected_plugins = Q2_ALLOWED_PLUGINS.copy()

for q2plugin, m in methods_to_add:
    inputs = m.signature.inputs.copy()
    outputs = m.signature.outputs.copy()
    parameters = m.signature.parameters.copy()
    add_method = True

    qname = q2plugin.name
    mid = m.id

    if qname in q2_expected_plugins:
        q2_expected_plugins.remove(qname)

    # storing this information in req_params so we can use internally
    # while calling call_qiime2
    req_params = {'qp-hide-plugin': ('string', qname),
                  'qp-hide-method': ('string', mid)}
    outputs_params = {}
    opt_params = {}
    to_delete = []
    for pname, element in inputs.items():
        qt_name, predicate = get_qiime2_type_name_and_predicate(element)

        if qt_name not in Q2_QIITA_SEMANTIC_TYPE:
            # As of qiime2-2022.2 this filters out:
            # Hierarchy
            #   gneiss dendrogram_heatmap
            #   gneiss ilr_hierarchical
            # SampleEstimator
            #   sample-classifier predict_classification
            #   sample-classifier predict_regression
            # ProcrustesStatistics
            #   emperor procrustes_plot
            add_method = False
            break

        add_method = True
        etype = Q2_QIITA_SEMANTIC_TYPE[qt_name]
        # these are special types as we can retrive internally
        if etype in 'phylogeny':
            ename = 'Phylogenetic tree'
            req_params[ename] = (
                'choice:["None", "Artifact tree, if exists"]', 'None')
            req_params['qp-hide-param' + ename] = ('string', pname)
            to_delete.append(pname)
        elif etype == 'FeatureData':
            if predicate == ['Taxonomy']:
                etype = 'FeatureData[Taxonomy]'
                ename = 'qp-hide-%s' % etype
                req_params[ename] = ('string', etype)
                to_delete.append(pname)
            elif predicate == ['Sequence']:
                req_params[element.description] = ('artifact', ['BIOM'])
                to_delete.append(pname)
            else:
                # As of qiime2-2022.2 this filters out:
                # predicate: Importance
                #   importance heatmap
                #   importances plot_feature_volatility
                # predicate: Differential
                #   differential ilr_phylogenetic_differential
                # predicate: AlignedSequence | Sequence
                #   data tabulate_seqs
                add_method = False
        elif etype == 'TaxonomicClassifier':
            default = qp_qiime2_dbs[0]
            qp_qiime2_dbs = ', '.join('"%s"' % db for db in qp_qiime2_dbs)

            ename = '%s (%s)' % (element.description, pname)
            req_params[ename] = ('choice:[%s]' % qp_qiime2_dbs, default)
            req_params['qp-hide-param' + ename] = ('string', pname)
        else:
            ename = f'{element.description} [{pname}]'

            # this is an odd one, first encountered:
            # feature-classifier fit-classifier-naive-bayes
            if ename == element.NOVALUE:
                # As of qiime2-2022.2 nothing is filtered here, so let's raise
                # an error so we can catch and review if this happens in the
                # future
                raise ValueError(f"[REVIEW] {pname} {mid} due to {ename}")
            else:
                req_params[ename] = ('artifact', [etype])
                req_params['qp-hide-param' + ename] = ('string', pname)

    if add_method:
        for td in to_delete:
            del inputs[td]

        # skipping if we know already that we will not add this method
        for pname, element in outputs.items():
            qt_name, predicate = get_qiime2_type_name_and_predicate(element)
            if (qt_name not in Q2_QIITA_SEMANTIC_TYPE or
                    qt_name in NOT_VALID_OUTPUTS):
                # As of qiime2-2022.2 this filters out:
                # Hierarchy
                #   gneiss assign_ids
                #   gneiss correlation_clustering
                #   gneiss gradient_clustering
                #   gneiss ilr_phylogenetic
                # Phylogeny
                #   gneiss ilr_phylogenetic_differential
                #   gneiss ilr_phylogenetic_ordination
                #   phylogeny align_to_tree_mafft_fasttree
                #   phylogeny align_to_tree_mafft_iqtree
                #   phylogeny align_to_tree_mafft_raxml
                #   phylogeny filter_tree
                # SampleEstimator
                #   longitudinal feature_volatility
                #   longitudinal maturity_index
                #   sample-classifier classify_samples
                #   sample-classifier fit_classifier
                #   sample-classifier fit_regressor
                #   sample-classifier regress_samples
                # ProcrustesStatistics
                #   diversity procrustes_analysis
                add_method = False
                break
            else:
                etype = Q2_QIITA_SEMANTIC_TYPE[qt_name]
                outputs_params[pname] = etype

        # we need to add the extra Qiita resulting table from assigning
        # taxonomy via classify_sklearn
        if qname == 'feature-classifier' and mid == 'classify_sklearn':
            outputs_params['Feature Table with Classification'] = 'BIOM'

    if not add_method:
        # for details of which commands are filtered due to this flag check
        # the code above
        continue

    total_inputs = len(inputs)
    if total_inputs == 0:
        # As of qiime2-2022.2 this filters out:
        # filtered_sequences filter_seqs
        continue

    for pname, element in parameters.items():
        tqt, predicate = get_qiime2_type_name_and_predicate(element)
        # there is a new primitive and we should raise an error
        if tqt not in PRIMITIVE_TYPES:
            raise ValueError(
                'There is a new type: %s, in %s %s (%s)' % (
                    tqt, qname, mid, pname))

        # predicate are the options for each parameter, note that it
        # can be a Choice/List or a Range (for Int/Floats). We ignore
        # numeric because they are ranges and we don't support ranges
        # in Qiita.
        # Note, the correct way to retrieve the latter is:
        # p.start, p.end, p.inclusive_start, p.inclusive_end
        # but for simplicity we will only retrive the
        data_type = PRIMITIVE_TYPES[tqt]
        default = element.default
        if (predicate is not None and PRIMITIVE_TYPES[tqt] not in (
                                      'float', 'integer')):
            vals = predicate.to_ast()['choices']
            data_type = 'choice:%s' % dumps(vals)
            default = vals[0]

        # if we are in the diversity plugin, the method starts with
        # alpha/beta, and the parameter is 'metric', we might want
        # to replace the technical names for user friendly ones
        value_pair = (mid, pname)
        if qname == 'diversity' and value_pair in RENAME_COMMANDS:
            # converting to list to the serialize doesn't complaint
            vals = list(RENAME_COMMANDS[value_pair])
            data_type = 'choice:%s' % dumps(vals)
            default = vals[0]

        # the diversity methods can have a choice param with no values
        # so we need to fix so users can actually select things;
        # however,we want to make sure that this is the only one,
        # if not, raise an error
        if data_type == 'choice' and default is None:
            # these are valid optional parameters
            valid_none_params = [
                ('emperor', 'plot'),
                ('gneiss', 'ilr_phylogenetic_ordination'),
                ('composition', 'ancombc')]
            if (qname, mid) in valid_none_params:
                data_type = 'string'
                default = ''
            else:
                error_msg = (
                    "There is an unexpected method (%s %s) with a "
                    "choice parameter (%s: %s), without default" % (
                        qname, mid, pname, element.description))
                raise ValueError(error_msg)

        if tqt == 'Metadata':
            # for filter_features we need to list the available filter qza
            if qname == 'feature-table' and mid == 'filter_features':
                qfq = ', '.join('"%s"' % qza for qza in qp_filtering_qza)
                ename = '%s (%s)' % (element.description, pname)
                opt_params[ename] = ('choice:["", %s]' % qfq, '')
                opt_params['qp-hide-metadata'] = ('string', ename)
            else:
                opt_params['qp-hide-metadata'] = ('string', pname)
        elif tqt == 'MetadataColumn':
            name = "Metadata column to use"
            opt_params[name] = ('string', '')
            name = 'qp-hide-param' + name
            opt_params[name] = ('string', 'qp-hide-metadata-field')
        else:
            ename = '%s (%s)' % (element.description, pname)

            if element.has_default():
                opt_params[ename] = (data_type, default)
                # we need to add the actual name of the parameter so we
                # can retrieve later
                opt_params['qp-hide-param' + ename] = ('string', pname)
            else:
                default = (default if default is not element.NOVALUE
                           else 'None')
                opt_params[ename] = (data_type, default)
                # we need to add the actual name of the parameter so we
                # can retrieve later
                opt_params['qp-hide-param' + ename] = ('string', pname)

    qiime_cmd = QiitaCommand("%s [%s]" % (m.name, mid), m.description,
                             call_qiime2, req_params, opt_params,
                             outputs_params, {'Default': {}},
                             analysis_only=True)

    plugin.register_command(qiime_cmd)

plugin.register_command(qiime_cmd)

if q2_expected_plugins:
    raise ValueError(f'Never saw plugin(s): {q2_expected_plugins}')
