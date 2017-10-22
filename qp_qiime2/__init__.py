# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------


from qiita_client import QiitaPlugin, QiitaCommand

from .qiime2 import (rarefy, beta_diversity, pcoa, beta_correlation,
                     alpha_diversity, alpha_correlation, taxa_barplot,
                     filter_samples, emperor, beta_group_significance)
from qiime2 import __version__ as qiime2_version


# Initialize the plugin
plugin = QiitaPlugin(
    'qiime2', qiime2_version, 'QIIME 2')

# Define the rarefy command
req_params = {
    'i-table': ('artifact', ['BIOM']),
    'p-sampling-depth': ['integer', 1000]
}
opt_params = {}
outputs = {'o-table': 'BIOM'}
dflt_param_set = {
    'Defaults': {
        'p-sampling-depth': 1000}
}
qiime_cmd = QiitaCommand(
    "Rarefy features", "Rarefy",
    rarefy, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define the beta_diversity command
req_params = {'i-table': ('artifact', ['BIOM'])}
opt_params = {
    'p-metric': [
        ('choice:["sokalsneath", "kulsinski", "mahalanobis", "matching", '
         '"euclidean", "correlation", "yule", "russellrao", "hamming", '
         '"jaccard", "braycurtis", "dice", "rogerstanimoto", "sqeuclidean", '
         '"cityblock", "sokalmichener", "cosine", "wminkowski", "seuclidean", '
         '"chebyshev", "canberra", "unweighted UniFrac", '
         '"weighted normalized UniFrac", "weighted unnormalized UniFrac"]'),
        'jaccard'],
    'i-tree': ['choice:["default", "None"]', 'None']}
outputs = {'distance_matrix': 'distance_matrix'}
dflt_param_set = {
    'Defaults': {
        'p-metric': 'jaccard',
        'i-tree': 'None'}
}
qiime_cmd = QiitaCommand(
    "Calculate beta diversity", "Beta Diversity",
    beta_diversity, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define the pcoa command
req_params = {'i-distance-matrix': ('artifact', ['distance_matrix'])}
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
req_params = {'i-distance-matrix': ('artifact', ['distance_matrix']),
              'm-metadata-category': ('string', '')}
opt_params = {'p-method': ['choice:["spearman", "pearson"]', 'spearman'],
              'p-permutations': ('integer', 999)}
outputs = {'q2_visualization': 'q2_visualization'}
dflt_param_set = {
    'Defaults': {
        'p-method': 'spearman',
        'p-permutations': 999}
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
        ('choice:["simpson", "enspie", "doubles", "mcintosh_d", "chao1", '
         '"kempton_taylor_q", "observed_otus", "singles", "strong", '
         '"robbins", "menhinick", "osd", "shannon", "brillouin_d", '
         '"margalef", "lladser_ci", "pielou_e", "michaelis_menten_fit", '
         '"ace", "chao1_ci", "goods_coverage", "esty_ci", "fisher_alpha", '
         '"berger_parker_d", "heip_e", "simpson_e", "dominance", '
         '"mcintosh_e", "lladser_pe", "gini_index", "faith_pd"]'),
        'observed_otus'],
    'Phylogenetic tree': ['choice:["default", "None"]', 'None']}
outputs = {'o-alpha-diversity': 'alpha_vector'}
dflt_param_set = {
    'Defaults': {
        'Diversity metric': 'observed_otus',
        'Phylogenetic tree': 'None'}
}
qiime_cmd = QiitaCommand(
    "Calculate alpha diversity", "Alpha Diversity",
    alpha_diversity, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define the alpha_correlation command
req_params = {'i-alpha-diversity': ('artifact', ['alpha_vector'])}
opt_params = {'p-method': ['choice:["spearman", "pearson"]', 'spearman']}
outputs = {'q2_visualization': 'q2_visualization'}
dflt_param_set = {
    'Defaults': {
        'p-method': 'spearman'}
}
qiime_cmd = QiitaCommand(
    "Calculate alpha correlation", "Alpha Corrrelation",
    alpha_correlation, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define the taxa barplot command
req_params = {'i-table': ('artifact', ['BIOM'])}
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
req_params = {'i-table': ('artifact', ['BIOM'])}
opt_params = {
    'p-min-frequency': ('integer', 0),
    'p-max-frequency': ('integer', 9223372036854775807),
    'p-min-features': ('integer', 0),
    'p-max-features': ('integer', 9223372036854775807),
    'p-where': ('string', '')}
outputs = {'o-table': 'BIOM'}
dflt_param_set = {
    'Defaults': {
        'p-min-frequency': 0,
        'p-max-frequency': 9223372036854775807,
        'p-min-features': 0,
        'p-max-features': 9223372036854775807,
        'p-where': ''}}
qiime_cmd = QiitaCommand(
    "Filter samples by metadata", "Filter Samples",
    filter_samples, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define the emperor command
req_params = {'i-pcoa': ('artifact', ['ordination_results'])}
opt_params = {'p-custom-axis': ('string', '')}
outputs = {'q2_visualization': 'q2_visualization'}
dflt_param_set = {'Defaults': {'p-custom-axis': ''}}
qiime_cmd = QiitaCommand(
    "Custom-axis Emperor plot", "Emperor plot",
    emperor, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)

# Define beta-group-significance command
req_params = {'i-distance-matrix': ('artifact', ['distance_matrix']),
              'm-metadata-category': ('string', '')}
opt_params = {'p-method': ['choice:["permanova", "anosim"]', 'permanova'],
              'p-pairwise': ['choice:["p-pairwise", "p-no-pairwise"]',
                             'p-pairwise'],
              'p-permutations': ('integer', 999)}
outputs = {'q2_visualization': 'q2_visualization'}
dflt_param_set = {
    'Defaults': {
        'p-method': 'permanova',
        'p-pairwise': 'p-pairwise',
        'p-permutations': 999}
}
qiime_cmd = QiitaCommand(
    "Calculate beta group significance", "Beta Group Significance",
    beta_group_significance, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)
