# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin, QiitaCommand

from .qiime2 import rarefy

__all__ = ['qiime2']


# Initialize the plugin
plugin = QiitaPlugin(
    'qiime2', '4.2017', 'QIIME 2')

# Define the rarefy command
req_params = {'i-table': ('artifact', ['BIOM'])}
opt_params = {
    'p-sampling-depth': ['integer', '10000']
}
outputs = {'o-table': 'BIOM'}
dflt_param_set = {
    'Defaults': {
        'p-sampling-depth': 1000}
}
qiime_cmd = QiitaCommand(
    "Rarefy", "Rarefy",
    rarefy, req_params, opt_params, outputs, dflt_param_set)
plugin.register_command(qiime_cmd)
