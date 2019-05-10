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
    PRIMITIVE_TYPES, call_qiime2, RENAME_COMMANDS)
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

# PLEASE READ:
# There are 2 main steps:
# 1. We are going loop over QIITA_Q2_SEMANTIC_TYPE (lookup table) so we can
# retrieve the q2plugin and their methods that work with that given Q2/Qiita
# semantic type. Then we will ignore any plugin not in Q2_ALLOWED_PLUGINS so
# we avoid adding plugins that we don't want; like deblur or dada2. Finally,
# we are going to loop over the different inputs, outputs and parameters from
# Q2 and convert them to QIITA's req_params, opt_params and outputs.
# 2. We are going to "force" add commands that are not added in the previous
# loop based on user demand
#
# Note that Qiita users like to have descriptions of the paramters
# (q2-description) vs. the parameter itself (q2-parameter) so to allow this
# we are going to store each parameter twice: one in the
# opt_params[q2-description]: value; and
# req_params['qp-hide-param' + q2-description]: q2-parameter

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
            # This currently filters (which are processing commands):
            # feature-classifier
            # fragment-insertion
            # quality-control
            # vsearch
            continue

        for m in methods:
            # after review of qiime2-2019.4 we decided to not add these methods
            if (q2plugin.name, m.id) in [('feature-table', 'group'),
                                         ('feature-table', 'filter_seqs')]:
                continue

            inputs = m.signature.inputs.copy()
            outputs = m.signature.outputs.copy()
            parameters = m.signature.parameters.copy()
            add_method = True

            # storing this information in req_params so we can use internally
            # while calling call_qiime2
            req_params = {'qp-hide-plugin': ('string', q2plugin.name),
                          'qp-hide-method': ('string', m.id)}
            outputs_params = {}
            opt_params = {}
            to_delete = []
            for pname, element in inputs.items():
                qt_name = element.qiime_type.to_ast()['name']
                if qt_name not in Q2_QIITA_SEMANTIC_TYPE:
                    add_method = False
                    break
                etype = Q2_QIITA_SEMANTIC_TYPE[qt_name]
                # these are special types as we can retrive internally
                if etype == 'phylogeny':
                    ename = 'Phylogenetic tree'
                    req_params[ename] = (
                        'choice:["None", "Artifact tree, if exists"]', 'None')
                    to_delete.append(pname)
                elif etype == 'FeatureData[Taxonomy]':
                    ename = 'qp-hide-%s' % etype
                    req_params[ename] = ('string', etype)
                    to_delete.append(pname)
                    # we are going to continue so we don't add this element
                    # twice
                    continue
                else:
                    ename = element.description
                    req_params[ename] = ('artifact', [etype])
                # we need to add the actual name of the parameter so we
                # can retrieve later
                req_params['qp-hide-param' + ename] = ('string', pname)

            for td in to_delete:
                del inputs[td]

            for pname, element in outputs.items():
                qt_name = element.qiime_type.to_ast()['name']
                if qt_name not in Q2_QIITA_SEMANTIC_TYPE:
                    add_method = False
                    break
                else:
                    etype = Q2_QIITA_SEMANTIC_TYPE[qt_name]
                    outputs_params[pname] = etype

            if len(inputs) != 1 or not add_method:
                # This is currently filtering out:
                # diversity mantel
                # diversity pcoa_biplot
                # diversity pcoa_biplot
                # diversity procrustes_analysis
                # emperor procrustes_plot
                # gneiss assign_ids
                # gneiss assign_ids
                # gneiss balance_taxonomy
                # gneiss balance_taxonomy
                # gneiss correlation_clustering
                # gneiss dendrogram_heatmap
                # gneiss gradient_clustering
                # gneiss gradient_clustering
                # gneiss ilr_hierarchical
                # gneiss ilr_phylogenetic
                # gneiss ilr_phylogenetic
                # longitudinal feature_volatility
                # longitudinal maturity_index
                # sample-classifier classify_samples
                # sample-classifier fit_classifier
                # sample-classifier fit_regressor
                # sample-classifier predict_classification
                # sample-classifier predict_regression
                # sample-classifier regress_samples
                # taxa filter_seqs
                continue

            for pname, element in parameters.items():
                tqt = element.qiime_type.to_ast()['name']
                # there is a new primitive and we should raise an error
                if tqt not in PRIMITIVE_TYPES:
                    raise ValueError(
                        'There is a new type: %s, in %s %s (%s)' % (
                            tqt, q2plugin.name, m.id, pname))

                # predicate are the options for each parameter, note that it
                # can be a Choice/List or a Range (for Int/Floats). We ignore
                # numeric because they are ranges and we don't support ranges
                # in Qiita.
                # Note, the correct way to retrieve the latter is:
                # p.start, p.end, p.inclusive_start, p.inclusive_end
                # but for simplicity we will only retrive the
                predicate = element.qiime_type.predicate
                data_type = PRIMITIVE_TYPES[tqt]
                default = element.default
                if (predicate is not None and PRIMITIVE_TYPES[tqt] not in (
                                              'float', 'integer')):
                    vals = predicate.to_ast()['choices']
                    data_type = 'choice:%s' % dumps(vals)
                    default = vals[0]

                qname = q2plugin.name
                mid = m.id
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
                    if qname == 'emperor' and mid == 'plot':
                        data_type = 'string'
                        default = ''
                    else:
                        error_msg = (
                            "There is an unexpected method (%s %s) with a "
                            "choice parameter (%s: %s), without default" % (
                                qname, mid, pname, element.description))
                        raise ValueError(error_msg)

                if tqt == 'Metadata':
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

            qiime_cmd = QiitaCommand("%s [%s]" % (m.name, m.id), m.description,
                                     call_qiime2, req_params, opt_params,
                                     outputs_params, {'Default': {}},
                                     analysis_only=True)

            plugin.register_command(qiime_cmd)

# 2. force adding extra commands
pm = PluginManager()

# 2.1 adding assign taxonomy
q2plugin = pm.plugins['feature-classifier']
m = q2plugin.methods['classify_sklearn']
qname = q2plugin.name
mid = m.id

req_params = {'qp-hide-plugin': ('string', q2plugin.name),
              'qp-hide-method': ('string', m.id)}
outputs_params = {}
opt_params = {}
for pname, element in m.signature.inputs.items():
    eqt = str(element.qiime_type)
    if eqt == 'FeatureData[Sequence]':
        req_params[element.description] = ('artifact', ['BIOM'])
    elif eqt == 'TaxonomicClassifier':
        default = qp_qiime2_dbs[0]
        qp_qiime2_dbs = ', '.join('"%s"' % db for db in qp_qiime2_dbs)

        ename = '%s (%s)' % (element.description, pname)
        req_params[ename] = ('choice:[%s]' % qp_qiime2_dbs, default)
        req_params['qp-hide-param' + ename] = ('string', pname)
    else:
        raise ValueError('Found unexpected input: "%s", in '
                         'feature-classifier classify_sklearn' % eqt)
for pname, element in m.signature.outputs.items():
    eqt = str(element.qiime_type)
    if eqt == 'FeatureData[Taxonomy]':
        outputs_params[pname] = eqt
        outputs_params['Feature Table with Classification'] = 'BIOM'
    else:
        raise ValueError('Found non expected output: "%s", in '
                         'feature-classifier classify_sklearn' % eqt)
for pname, element in m.signature.parameters.items():
    tqt = element.qiime_type.to_ast()['name']
    if tqt not in PRIMITIVE_TYPES:
        raise ValueError(
            'There is a new type: %s, in %s %s (%s)' % (
                tqt, q2plugin.name, m.id, pname))
    predicate = element.qiime_type.predicate
    data_type = PRIMITIVE_TYPES[tqt]
    default = element.default
    if (predicate is not None and PRIMITIVE_TYPES[tqt] not in (
                                  'float', 'integer')):
        vals = predicate.to_ast()['choices']
        data_type = 'choice:%s' % dumps(vals)
        default = vals[0]

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

qiime_cmd = QiitaCommand("%s [%s]" % (m.name, m.id), m.description,
                         call_qiime2, req_params, opt_params,
                         outputs_params, {'Default': {}}, analysis_only=True)
plugin.register_command(qiime_cmd)
