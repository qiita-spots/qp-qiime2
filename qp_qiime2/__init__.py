# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------


from qiita_client import QiitaPlugin, QiitaCommand

from .qiime2 import rarefy, beta_diversity
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
outputs = {'distance-matrix': 'distance_matrix'}
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
