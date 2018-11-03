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
    QIITA_Q2_SEMANTIC_TYPE, Q2_QIITA_SEMANTIC_TYPE, Q2_ALLOWED_PLUGINS,
    PRIMITIVE_TYPES, call_qiime2, ALPHA_DIVERSITY_METRICS,
    ALPHA_DIVERSITY_METRICS_PHYLOGENETIC, BETA_DIVERSITY_METRICS,
    BETA_DIVERSITY_METRICS_PHYLOGENETIC,
    BETA_DIVERSITY_METRICS_PHYLOGENETIC_ALT)
from qiime2 import __version__ as qiime2_version
from qiime2.sdk.util import actions_by_input_type


# Initialize the qiita_plugin
plugin = QiitaPlugin('qiime2', qiime2_version, 'QIIME 2')

# PLEASE READ:
# We are going loop over QIITA_Q2_SEMANTIC_TYPE (lookup table)
# so we can retrieve the q2plugin and their methods that work with that
# given Q2/Qiita semantic type. Then we will ignore any plugin not in
# Q2_ALLOWED_PLUGINS so we avoid adding plugins that we don't want; like
# deblur or dada2. Finally, we are going to loop over the different inputs,
# outputs and parameters from Q2 and convert them to QIITA's req_params,
# opt_params and outputs.
# Also, note that Qiita users like to have descriptions of the paramters
# (q2-description) vs. the parameter itself (q2-parameter) so to allow this
# we are going to store each parameter twice: one in the
# opt_params[q2-description]: value; and
# req_params['qp-hide-param' + q2-description]: q2-parameter
for qiita_artifact, q2_artifact in QIITA_Q2_SEMANTIC_TYPE.items():
    for q2plugin, methods in actions_by_input_type(str(q2_artifact)):
        # note that the qiita_artifact are strings not objects
        if qiita_artifact.startswith('BIOM'):
            qiita_artifact = 'BIOM'

        if q2plugin.name not in Q2_ALLOWED_PLUGINS:
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
                if element.qiime_type not in Q2_QIITA_SEMANTIC_TYPE:
                    add_method = False
                    break
                etype = Q2_QIITA_SEMANTIC_TYPE[element.qiime_type]
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
                    # we are gonna continue so we don't add this element twice
                    continue
                else:
                    ename = element.description
                    req_params[ename] = ('artifact', [etype])
                # we need to add the actual name of the parameter so we
                # can retrieve later
                req_params['qp-hide-param' + ename] = ('string', pname)

            outputs_params = {}
            for pname, element in outputs.items():
                if element.qiime_type not in Q2_QIITA_SEMANTIC_TYPE:
                    add_method = False
                    break
                else:
                    etype = Q2_QIITA_SEMANTIC_TYPE[element.qiime_type]
                    if etype.startswith('BIOM'):
                        etype = 'BIOM'
                    # this one is to "fix" the templates for phylogenetic
                    # methods, like phylogenetic_distance_matrix
                    elif etype.startswith('phylogenetic_'):
                        etype = etype[len('phylogenetic_'):]
                    edesc = (element.description if element.has_description()
                             else pname)
                    outputs_params[edesc] = etype

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
                    vals = list(predicate.iter_boundaries())
                    data_type = 'choice:%s' % dumps(vals)
                    default = vals[0]
                # the diversity methods can have a choice param with no values
                # so we need to fix so users can actually select things;
                # however,we want to make sure that this is the only one,
                # if not, raise an error
                if data_type == 'choice' and default is None:
                    qname = q2plugin.name
                    mid = m.id
                    error_msg = ("There is an unexpected method (%s %s) with "
                                 "a choice parameter (%s), without default" % (
                                    qname, mid, element.description))
                    # if we are in the diversity choice option, we might want
                    # to replace the "nerd" names for user friendly ones
                    if qname == 'diversity':
                        am = set(ALPHA_DIVERSITY_METRICS)
                        bm = set(BETA_DIVERSITY_METRICS)
                        amp = set(ALPHA_DIVERSITY_METRICS_PHYLOGENETIC)
                        bmp = set(BETA_DIVERSITY_METRICS_PHYLOGENETIC)
                        bmpa = set(BETA_DIVERSITY_METRICS_PHYLOGENETIC_ALT)
                        vals = {
                            'alpha': am,
                            'alpha_phylogenetic': amp,
                            'beta': bm,
                            'beta_phylogenetic': bmp,
                            'beta_phylogenetic_alt': bmpa,
                            'alpha_rarefaction': am.union(amp),
                            'beta_rarefaction': bm.union(bmp)
                        }
                        if mid not in mid:
                            raise ValueError(error_msg)
                        # converting to list to the serialize doesn't complaint
                        vals = list(vals[mid])
                        data_type = 'choice:%s' % dumps(vals)
                        default = vals[0]
                    elif qname == 'emperor' and mid == 'plot':
                        data_type = 'string'
                        default = ''
                    else:
                        raise ValueError(error_msg)

                if pname == 'metadata':
                    # Q2 does some CLI magic when dealing with mapping
                    # files, if the method requires a column, the CLI will
                    # request the filepath, parse it and then pass as metadata
                    # column. For fun, both cases full/column metadata are
                    # called metadata. However, for Qiita we will need to make
                    # this difference more obvious.
                    if data_type == 'string':
                        # as this one needs input from the user, we will create
                        # as any other opt_params
                        name = "Metadata column to use"
                        opt_params[name] = ('string', '')
                        name = 'qp-hide-param' + name
                        opt_params[name] = ('string', 'qp-hide-metadata-field')
                    else:
                        opt_params['qp-hide-metadata'] = ('string', pname)
                else:
                    ename = element.description
                    if element.has_default():
                        opt_params[ename] = (data_type, default)
                        # we need to add the actual name of the parameter so we
                        # can retrieve later
                        opt_params['qp-hide-param' + ename] = ('string', pname)
                    else:
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
