# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

import traceback
import sys
from os.path import join, dirname, abspath
from os import environ, close
from configparser import ConfigParser
from tempfile import mkstemp

from qiita_client import QiitaClient, ArtifactInfo
from q2_types import (
    FeatureTable, Frequency, RelativeFrequency, PresenceAbsence)

from biom import load_table, Table
from biom.util import biom_open
from qiime.sdk import PluginManager, Artifact


def get_plugin_method_lookup():
    plugin_manager = PluginManager()
    result = {}
    for plugin_name, plugin in plugin_manager.plugins.items():
        for method_name, method in plugin.methods.items():
            key = '%s: %s' % (plugin_name, method_name)
            result[key] = method
    return result


def biom_artifact_input_translator(filepaths):
    biom_fp = filepaths['biom'][0]
    return load_table(biom_fp)


def biom_artifact_output_translator(artifact):
    biom_table = artifact.view(Table)
    fd, temp_file_name = mkstemp(suffix=".biom")
    close(fd)
    with biom_open(temp_file_name, 'w') as f:
        biom_table.to_hdf5(f, "QIITA-QIIME 2 plugin")
    return temp_file_name, 'biom'


def get_artifact_translators_lookup():
    # TODO (GREG): Figure out how to lookup all subtypes w/o explicit definitions
    # TODO (GREG): make FeatureTable, etc, hashable see issue qiime2 #46
    biom_tuple = (biom_artifact_input_translator,
                  biom_artifact_output_translator,
                  'BIOM')
    return {FeatureTable: biom_tuple,
            FeatureTable[Frequency]: biom_tuple,
            FeatureTable[RelativeFrequency]: biom_tuple,
            FeatureTable[PresenceAbsence]: biom_tuple}


def retrieve_artifact_info(qclient, artifact_id, artifact_type):
    a_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    filepaths = a_info['files']
    translator_lookup = get_artifact_translators_lookup()
    input_translator, _, _ = translator_lookup[artifact_type[0]]
    artifact_data = input_translator(filepaths)
    # TODO: None is the provenance
    a = Artifact._from_view(artifact_data, artifact_type[0], None)
    return a


def job_id_to_method(qclient, job_id, job_info):
    method_lookup = get_plugin_method_lookup()
    method = method_lookup[job_info['command']]
    return method


def format_artifacts_info(output_names, output_artifacts):
    result = []
    translator_lookup = get_artifact_translators_lookup()

    for output_name, a in zip(output_names, output_artifacts):
        _, output_translator, qiita_type = translator_lookup[a.type]
        out_file, fp_type = output_translator(a)
        result.append(
            ArtifactInfo(output_name, qiita_type, [(out_file, fp_type)]))
        # TODO: figure out how to hadle the case where multiple qiime artifacts
        # map to one Qiita artifact

    return result


def execute_task(qclient, job_id, job_info):
    method = job_id_to_method(qclient, job_id, job_info)

    parameters = job_info['parameters']

    # execute workflow
    inputs = {
        ia_name: retrieve_artifact_info(qclient, parameters[ia_name], ia_type)
        for ia_name, ia_type in method.signature.inputs.items()}
    input_parameters = {}
    input_parameters = {ip_name: parameters[ip_name]
                        for ip_name in method.signature.parameters}

    args = inputs
    args.update(input_parameters)
    # Method can return either None, a single artifact or a tuple of artifacts
    output_artifacts = method(**args)

    if output_artifacts is None:
        return None
    elif isinstance(output_artifacts, Artifact):
        output_artifacts = (output_artifacts, )

    artifacts_info = format_artifacts_info(method.signature.outputs,
                                           output_artifacts)

    return True, artifacts_info, ""


def execute_job(url, job_id, output_dir):
    # Set up the QiitaClient - Probably another function?
    dflt_conf_fp = join(dirname(abspath(__file__)), 'support_files',
                        'config_file.cfg')
    conf_fp = environ.get('QP_QIIME2_CONFIG_FP', dflt_conf_fp)
    config = ConfigParser()
    with open(conf_fp, 'U') as conf_file:
        config.readfp(conf_file)

    qclient = QiitaClient(url, config.get('main', 'CLIENT_ID'),
                          config.get('main', 'CLIENT_SECRET'),
                          server_cert=config.get('main', 'SERVER_CERT'))

    # Request job information. If there is a problem retrieving the job
    # information, the QiitaClient already raises an error
    job_info = qclient.get_job_info(job_id)
    # Starting the heartbeat
    qclient.start_heartbeat(job_id)

    try:
        success, artifacts_info, error_msg = execute_task(qclient, job_id,
                                                          job_info)
    except Exception:
        exc_str = ''.join(traceback.format_exception(*sys.exc_info()))
        error_msg = ("Error executing %s:\n%s" % (job_info['command'],
                                                  exc_str))
        success = False
        artifacts_info = None

    # The job completed
    qclient.complete_job(job_id, success, error_msg=error_msg,
                         artifacts_info=artifacts_info)

    print(error_msg)
