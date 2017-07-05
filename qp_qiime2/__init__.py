# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------


from qiita_client import QiitaPlugin, QiitaCommand

from .qiime2 import rarefy, beta_diversity, pcoa, beta_correlation
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
    "Rarefy", "Rarefy",
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
    "beta_diversity", "Beta Diversity",
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
    "pcoa", "Principal Coordinate Analysis",
    pcoa, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)


# Define the beta_correlation command
req_params = {'i-distance-matrix': ('artifact', ['distance_matrix']),
              'm-metadata-file': ('artifact', ['BIOM'])}
opt_params = {'m-metadata-category': ('string', ''),
              'p-method': ['choice:["spearman", "pearson"]', 'spearman'],
              'p-permutations': ('integer', 999)}
# outputs = {'qiime2-visualization': 'qiime2-visualization'}
outputs = {'q2_visualization': 'q2_visualization'}
dflt_param_set = {
    'Defaults': {
        'm-metadata-category': '#SampleID',
        'p-method': 'spearman',
        'p-permutations': 999}
}
qiime_cmd = QiitaCommand(
    "beta_correlation", "Beta Corrrelation",
    beta_correlation, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)
