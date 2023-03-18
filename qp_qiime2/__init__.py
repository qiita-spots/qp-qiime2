# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin

from qiime2 import __version__ as qiime2_version
from qiime2.sdk import PluginManager
from qiime2.sdk.util import actions_by_input_type

from .qp_qiime2 import (
    QIITA_Q2_SEMANTIC_TYPE, Q2_ANALYSIS_PLUGINS, Q2_PROCESSING_PLUGINS,
    Q2_EXTRA_COMMANDS)
from .util import register_qiime2_commands


# Initialize the qiita_plugin
plugin = QiitaPlugin('qiime2', qiime2_version, 'QIIME 2 - Analysis')


# PLEASE READ:
# There are 2 main steps:
# 1. We are going to loop on the Q2_EXTRA_COMMANDS to check which extra
# commands we want to add, extra from those given in the next step
# 2. We are going loop over QIITA_Q2_SEMANTIC_TYPE (lookup table) so we can
# retrieve the q2plugin and their methods that work with that given Q2/Qiita
# semantic type. Then we will ignore any plugin not in Q2_ANALYSIS_PLUGINS so
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

        if q2plugin.name not in Q2_ANALYSIS_PLUGINS:
            # As of qiime2-2022.11 this filters out:
            # alignment
            # deblur
            # diversity-lib
            # feature-classifier
            # fragment-insertion
            # greengenes2
            # quality-control
            # sourcetracker2
            # vsearch
            continue

        for m in methods:
            # after review of qiime2-2019.4 we decided to not add these methods
            if (q2plugin.name, m.id) not in [('feature-table', 'group'),
                                             ('feature-table', 'filter_seqs'),
                                             # qiime2-2022.11 we added this:
                                             ('composition', 'ancombc')]:
                methods_to_add.append((q2plugin, m))

# make sure we have seen all expected analysis plugins
q2_expected_plugins = register_qiime2_commands(
    plugin, methods_to_add, Q2_ANALYSIS_PLUGINS.copy())
if q2_expected_plugins:
    raise ValueError(f'Never saw plugin(s): {q2_expected_plugins}')

# make sure we have seen all expected processing plugins
gg2 = pm.plugins['greengenes2']
methods = [
    (gg2, gg2.methods['filter_features']),
    (gg2, gg2.actions['non_v4_16s']),
]
q2_expected_plugins = register_qiime2_commands(
    plugin, methods, Q2_PROCESSING_PLUGINS.copy(), False)
if q2_expected_plugins:
    raise ValueError(f'Never saw plugin(s): {q2_expected_plugins}')
