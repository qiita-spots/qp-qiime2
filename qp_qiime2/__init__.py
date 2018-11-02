# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from json import dumps

from qiita_client import QiitaPlugin, QiitaCommand

from .qp_qiime2 import (
    QIITA_Q2_ARTIFACTS, Q2_QIITA_ARTIFACTS, PLUGINS, PRIMITIVE_TYPES,
    call_qiime2, ALPHA_DIVERSITY_METRICS)
from qiime2 import __version__ as qiime2_version
from qiime2.sdk.util import actions_by_input_type


# Initialize the qiita_plugin
plugin = QiitaPlugin('qiime2', qiime2_version, 'QIIME 2')

all_tbla = []
for qiita_artifact, q2_artifact in QIITA_Q2_ARTIFACTS.items():
    for q2plugin, methods in actions_by_input_type(str(q2_artifact)):
        if qiita_artifact.startswith('BIOM'):
            qiita_artifact = 'BIOM'

        if q2plugin.name not in PLUGINS[qiita_artifact]:
            # This currently filters (which are processing commands):
            # feature-classifier
            # quality-control
            # vsearch
            # fragment-insertion
            continue
        for m in methods:
            inputs = m.signature.inputs.copy()
            outputs = m.signature.outputs.copy()
            parameters = m.signature.parameters.copy()
            add_method = True

            # storing this information in req_params so we can use internally
            # while calling call_qiime2
            req_params = {'qp-hide-plugin': ('string', q2plugin.name),
                          'qp-hide-method': ('string', m.id)}
            for pname, element in inputs.items():
                if element.qiime_type not in Q2_QIITA_ARTIFACTS:
                    add_method = False
                    break
                etype = Q2_QIITA_ARTIFACTS[element.qiime_type]
                if etype.startswith('BIOM'):
                    etype = 'BIOM'

                # these are special types as we can retrive internally
                if etype == 'phylogeny':
                    ename = 'Phylogenetic tree'
                    req_params[ename] = (
                        'choice:["default", "None"]', 'None')
                    # deleting so we don't count it as part of the inputs
                    del inputs[pname]
                elif etype == 'taxonomy':
                    ename = 'qp-hide-%s' % etype
                    req_params[ename] = ('string', etype)
                    # deleting so we don't count it as part of the inputs
                    del inputs[pname]
                else:
                    ename = element.description
                    req_params[ename] = ('artifact', [etype])
                # we need to add the actual name of the parameter so we
                # can retrieve later
                req_params['qp-hide-param' + ename] = ('string', pname)

            outputs_params = {}
            for pname, element in outputs.items():
                if element.qiime_type not in Q2_QIITA_ARTIFACTS:
                    add_method = False
                    break
                else:
                    etype = Q2_QIITA_ARTIFACTS[element.qiime_type]
                    if etype.startswith('BIOM'):
                        etype = 'BIOM'
                    elif etype.startswith('phylogenetic_'):
                        etype = etype[len('phylogenetic_'):]
                    edesc = (element.description if element.has_description()
                             else pname)
                    outputs_params[edesc] = etype
                    # Todo: how to save the name of the outputs??
                    # # we need to add the actual name of the parameter so we
                    # # can retrieve later
                    # req_params['qp-hide-param_' + ename] = ('string', pname)

            if len(inputs) != 1 or not add_method:
                # This is currently filtering out:
                # emperor procrustes_plot
                # diversity procrustes_analysis
                # diversity pcoa_biplot
                # diversity mantel
                # longitudinal first_distances
                # sample-classifier classify_samples_from_dist
                # diversity pcoa_biplot
                # gneiss gradient_clustering
                # feature-table summarize
                # feature-table presence_absence
                # longitudinal plot_feature_volatility
                # longitudinal first_differences
                # gneiss ilr_phylogenetic
                # gneiss correlation_clustering
                # gneiss assign_ids
                # gneiss gradient_clustering
                # gneiss dendrogram_heatmap
                # gneiss balance_taxonomy
                # gneiss ilr_hierarchical
                # feature-table filter_seqs
                # feature-table summarize
                # feature-table presence_absence
                # longitudinal maturity_index
                # longitudinal feature_volatility
                # phylogeny filter_table
                # sample-classifier classify_samples
                # sample-classifier predict_regression
                # sample-classifier fit_regressor
                # sample-classifier fit_classifier
                # sample-classifier regress_samples_ncv
                # sample-classifier regress_samples
                # sample-classifier classify_samples_ncv
                # sample-classifier predict_classification
                # composition add_pseudocount
                # feature-table summarize
                continue

            opt_params = {}
            for pname, element in parameters.items():
                tqt = type(element.qiime_type)
                # there is a new primitive and we should raise an error
                if tqt not in PRIMITIVE_TYPES:
                    raise ValueError(
                        'There is a new type: %s' % element.qiime_type)

                # predicate are the options for each parameter, note that it
                # can be a Choise/List or a Range (for Int/Floats). We ignore
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
                    vals = list(predicate.iter_boundaries())
                    data_type = 'choice:%s' % dumps(vals)
                    default = vals[0]
                # alpha_rarefaction can have a choice param with no values so
                # we need to fix so users can actually select things; however,
                # we want to make sure that this is the only one, if not, raise
                # an error
                if data_type == 'choice' and default is None:
                    qname = q2plugin.name
                    mid = m.id
                    if qname == 'diversity' and mid == 'alpha_rarefaction':
                        vals = list(ALPHA_DIVERSITY_METRICS.keys())
                        data_type = 'choice:%s' % dumps(vals)
                        default = vals[0]
                    elif qname == 'emperor' and mid == 'plot':
                        vals = list(ALPHA_DIVERSITY_METRICS.keys())
                        data_type = 'string'
                        default = 'None'
                    else:
                        raise ValueError(
                            "There is an unexpected method (%s %s) with a "
                            "choice parameter (%s), without default" % (
                                qname, mid, element.description))

                # only optional parameters have defaults
                if element.has_default():
                    if pname == 'metadata':
                        if data_type == 'string':
                            opt_params['qp-hide-metadata-field'] = (
                                'string', '')
                        elif data_type == 'mapping':
                            opt_params['qp-hide-metadata'] = ('string', '')
                        else:
                            raise ValueError(
                                "Not valid metadata data type: %s" % data_type)
                    else:
                        ename = element.description
                        opt_params[ename] = (data_type, default)
                        # we need to add the actual name of the parameter so we
                        # can retrieve later
                        opt_params['qp-hide-param' + ename] = ('string', pname)
                else:
                    if pname == 'metadata':
                        if data_type == 'string':
                            req_params['qp-hide-metadata-field'] = (
                                'string', '')
                        elif data_type == 'mapping':
                            req_params['qp-hide-metadata'] = ('string', '')
                        else:
                            raise ValueError(
                                "Not valid metadata data type: %s" % data_type)
                    else:
                        ename = element.description
                        default = (default if default is not element.NOVALUE
                                   else 'None')
                        req_params[ename] = (data_type, default)
                        # we need to add the actual name of the parameter so we
                        # can retrieve later
                        req_params['qp-hide-param' + ename] = ('string', pname)

            qiime_cmd = QiitaCommand(
                m.name, m.description, call_qiime2, req_params, opt_params,
                outputs_params, {'Defaut': {}}, analysis_only=True)

            plugin.register_command(qiime_cmd)
