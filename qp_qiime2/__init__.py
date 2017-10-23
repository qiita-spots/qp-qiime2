# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from json import dumps

from qiita_client import QiitaPlugin, QiitaCommand

from .qiime2 import (rarefy, beta_diversity, pcoa, beta_correlation,
                     alpha_diversity, alpha_correlation, taxa_barplot,
                     filter_samples, emperor, beta_group_significance,
                     ALPHA_DIVERSITY_METRICS, BETA_DIVERSITY_METRICS,
                     ALPHA_CORRELATION_METHODS, BETA_CORRELATION_METHODS,
                     BETA_GROUP_SIG_METHODS, BETA_GROUP_SIG_TYPE)
from qiime2 import __version__ as qiime2_version


# Initialize the plugin
plugin = QiitaPlugin(
    'qiime2', qiime2_version, 'QIIME 2')

# Define the rarefy command
req_params = {
    'BIOM table': ('artifact', ['BIOM']),
    'Sampling depth': ['integer', 1000]
}
opt_params = {}
outputs = {'o-table': 'BIOM'}
dflt_param_set = {
    'Defaults': {
        'Sampling depth': 1000}
}
qiime_cmd = QiitaCommand(
    "Rarefy features", "Rarefy",
    rarefy, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define the beta_diversity command
req_params = {'BIOM table': ('artifact', ['BIOM'])}
opt_params = {
    'Diversity metric': [
        'choice:%s' % dumps(list(BETA_DIVERSITY_METRICS.keys())),
        'Jaccard similarity index'],
    'Phylogenetic tree': ['choice:["default", "None"]', 'None']}
outputs = {'distance_matrix': 'distance_matrix'}
dflt_param_set = {
    'Defaults': {
        'Diversity metric': 'Jaccard similarity index',
        'Phylogenetic tree': 'None'}
}
qiime_cmd = QiitaCommand(
    "Calculate beta diversity", "Beta Diversity",
    beta_diversity, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define the pcoa command
req_params = {'Distance matrix': ('artifact', ['distance_matrix'])}
opt_params = {}
outputs = {'o-pcoa': 'ordination_results'}
dflt_param_set = {
    'Defaults': {}
}
qiime_cmd = QiitaCommand(
    "Generate principal coordinates analysis (PCoA)",
    "Principal Coordinate Analysis",
    pcoa, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define the beta_correlation command
req_params = {'Distance matrix': ('artifact', ['distance_matrix']),
              'Metadata category': ('string', '')}
opt_params = {'Correlation method':
              ['choice:%s' % dumps(list(BETA_CORRELATION_METHODS)),
               'Spearman'],
              'Number of permutations': ('integer', 999)}
outputs = {'q2_visualization': 'q2_visualization'}
dflt_param_set = {
    'Defaults': {
        'Correlation method': 'spearman',
        'Number of permutations': 999}
}
qiime_cmd = QiitaCommand(
    "Calculate beta correlation", "Beta Corrrelation",
    beta_correlation, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define the alpha command
req_params = {'BIOM table': ('artifact', ['BIOM'])}
opt_params = {
    'Diversity metric': [
        'choice:%s' % dumps(list(ALPHA_DIVERSITY_METRICS.keys())),
        'Number of distinct features'],
    'Phylogenetic tree': ['choice:["default", "None"]', 'None']}
outputs = {'o-alpha-diversity': 'alpha_vector'}
dflt_param_set = {
    'Defaults': {
        'Diversity metric': 'Number of distinct features',
        'Phylogenetic tree': 'None'}
}
qiime_cmd = QiitaCommand(
    "Calculate alpha diversity", "Alpha Diversity",
    alpha_diversity, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define the alpha_correlation command
req_params = {'Alpha vectors': ('artifact', ['alpha_vector'])}
opt_params = {'Correlation method':
              ['choice:%s' % dumps(list(ALPHA_CORRELATION_METHODS)),
               'Spearman']}
outputs = {'q2_visualization': 'q2_visualization'}
dflt_param_set = {'Defaults': {'Correlation method': 'Spearman'}}
qiime_cmd = QiitaCommand(
    "Calculate alpha correlation", "Alpha Corrrelation",
    alpha_correlation, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define the taxa barplot command
req_params = {'BIOM table': ('artifact', ['BIOM'])}
opt_params = {}
outputs = {'q2_visualization': 'q2_visualization'}
dflt_param_set = {
    'Defaults': {}
}
qiime_cmd = QiitaCommand(
    "Summarize taxa", "Taxa Barplot",
    taxa_barplot, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define the filtering samples from biom command
req_params = {'BIOM table': ('artifact', ['BIOM'])}
opt_params = {
    'Minimum feature frequency across samples': ('integer', 0),
    'Maximum feature frequency across samples':
        ('integer', 9223372036854775807),
    'Minimum features per sample': ('integer', 0),
    'Maximum features per sample': ('integer', 9223372036854775807),
    'SQLite WHERE-clause': ('string', '')}
outputs = {'o-table': 'BIOM'}
dflt_param_set = {
    'Defaults': {
        'Minimum feature frequency across samples': 0,
        'Maximum feature frequency across samples': 9223372036854775807,
        'Minimum features per sample': 0,
        'Maximum features per sample': 9223372036854775807,
        'SQLite WHERE-clause': ''}}
qiime_cmd = QiitaCommand(
    "Filter samples by metadata", "Filter Samples",
    filter_samples, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define the emperor command
req_params = {'Ordination results': ('artifact', ['ordination_results'])}
opt_params = {'Custom axis': ('string', '')}
outputs = {'q2_visualization': 'q2_visualization'}
dflt_param_set = {'Defaults': {'Custom axis': ''}}
qiime_cmd = QiitaCommand(
    "Custom-axis Emperor plot", "Emperor plot",
    emperor, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define beta-group-significance command
req_params = {'Distance matrix': ('artifact', ['distance_matrix']),
              'Metadata category': ('string', '')}
opt_params = {'Method':
              ['choice:%s' % dumps(list(BETA_GROUP_SIG_METHODS.keys())),
               'PERMANOVA'],
              'Comparison type':
              ['choice:%s' % dumps(list(BETA_GROUP_SIG_TYPE.keys())),
               'Pairwise'],
              'Number of permutations': ('integer', 999)}
outputs = {'q2_visualization': 'q2_visualization'}
dflt_param_set = {
    'Defaults': {
        'Method': 'permanova',
        'Comparison type': 'p-pairwise',
        'Number of permutations': 999}
}
qiime_cmd = QiitaCommand(
    "Calculate beta group significance", "Beta Group Significance",
    beta_group_significance, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)
