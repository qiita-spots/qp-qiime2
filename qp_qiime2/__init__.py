# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------


from qiita_client import QiitaPlugin, QiitaCommand

from .qiime2 import rarefy, beta_diversity, pcoa, alpha_diversity
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

# Define the alpha command
req_params = {'i-table': ('artifact', ['BIOM'])}
opt_params = {
    'p-metric': [
        ('choice:["simpson", "enspie", "doubles", "mcintosh_d", "chao1", '
         '"kempton_taylor_q", "observed_otus", "singles", "strong", '
         '"robbins", "menhinick", "osd", "shannon", "brillouin_d", '
         '"margalef", "lladser_ci", "pielou_e", "michaelis_menten_fit", '
         '"ace", "chao1_ci", "goods_coverage", "esty_ci", "fisher_alpha", '
         '"berger_parker_d", "heip_e", "simpson_e", "dominance", '
         '"mcintosh_e", "lladser_pe", "gini_index", "faith_pd"]'),
        'observed_otus'],
    'i-tree': ['choice:["default", "None"]', 'None']}
#
# change!
#
outputs = {'o-alpha-diversity': 'ordination_results'}
dflt_param_set = {
    'Defaults': {
        'p-metric': 'observed_otus',
        'i-tree': 'None'}
}
qiime_cmd = QiitaCommand(
    "alpha_diversity", "Alpha Diversity",
    alpha_diversity, req_params, opt_params, outputs, dflt_param_set,
    analysis_only=True)
plugin.register_command(qiime_cmd)
