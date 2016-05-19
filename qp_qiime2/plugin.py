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

from qiita_client import QiitaClient
from q2_types import (
    FeatureTable, Frequency, RelativeFrequency, PresenceAbsence)

from biom import load_table
from biom.util import biom_open
from qiime.sdk import PluginManager, Artifact, SubprocessExecutor


def get_plugin_workflow_lookup():
    plugin_manager = PluginManager()
    result = {}
    for plugin_name, plugin in plugin_manager.plugins.items():
        for workflow_name, workflow in plugin.workflows.items():
            key = '%s: %s' % (plugin_name, workflow_name)
            result[key] = workflow
    return result


def biom_artifact_input_translator(filepaths):
    # TODO (JOSE): Change filepaths to be a dictionary - ongoing work - qiita
    biom_fp = None

    for fp, fp_type in filepaths:
        if fp_type == 'biom':
            biom_fp = fp
            break

    return load_table(biom_fp)


def biom_artifact_output_translator(artifact):
    biom_table = artifact.data
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
    return {str(FeatureTable): biom_tuple,
            str(FeatureTable[Frequency]): biom_tuple,
            str(FeatureTable[RelativeFrequency]): biom_tuple,
            str(FeatureTable[PresenceAbsence]): biom_tuple}


def retrieve_artifact_info(qclient, artifact_id, artifact_type):
    fps_info = qclient.get("/qiita_db/artifacts/%s/filepaths/" % artifact_id)
    filepaths = fps_info['filepaths']
    translator_lookup = get_artifact_translators_lookup()
    input_translator, _, _ = translator_lookup[str(artifact_type)]
    artifact_data = input_translator(filepaths)
    fd, temp_file_name = mkstemp(suffix='.qtf')
    close(fd)
    # TODO: None is the provenance
    Artifact.save(artifact_data, artifact_type, None, temp_file_name)
    return temp_file_name


def job_id_to_workflow_executor_arguments(qclient, job_id, job_info):
    wf_lookup = get_plugin_workflow_lookup()
    wf = wf_lookup[job_info['command']]
    parameters = job_info['parameters']

    input_artifact_fps = {}
    for ia_name, ia_type in wf.signature.input_artifacts.items():
        artifact_fp = retrieve_artifact_info(qclient, parameters[ia_name],
                                             ia_type)
        input_artifact_fps[ia_name] = artifact_fp

    input_parameters = {ip_name: parameters[ip_name]
                        for ip_name in wf.signature.input_parameters}

    output_artifacts = {}
    for name in wf.signature.output_artifacts.keys():
        fd, temp_file_name = mkstemp(suffix='.qtf')
        close(fd)
        output_artifacts[name] = temp_file_name

    return wf, input_artifact_fps, input_parameters, output_artifacts


def format_artifacts_info(output_artifacts):
    result = []
    translator_lookup = get_artifact_translators_lookup()
    for output_name, fp in output_artifacts.items():
        a = Artifact(fp)
        _, output_translator, qiita_type = translator_lookup[str(a.type)]
        out_file, fp_type = output_translator(a)
        result.append([output_name, qiita_type, [(out_file, fp_type)]])
        # TODO: figure out how to hadle the case where multiple qiime artifacts
        # map to one Qiita artifact

    return result


def execute_task(qclient, job_id, job_info):
    wf, input_artifact_fps, input_parameters, output_artifacts = \
        job_id_to_workflow_executor_arguments(qclient, job_id, job_info)

    # execute workflow
    executor = SubprocessExecutor()
    future_ = executor(wf,
                       input_artifact_fps,
                       input_parameters,
                       output_artifacts)

    completed_process = future_.result()
    if completed_process.returncode == 0:
        success = True
        error_msg = ""
        artifacts_info = format_artifacts_info(output_artifacts)
    else:
        success = False
        error_msg = "StdOut: %s\nStdErr: %s" % (completed_process.stdout,
                                                completed_process.stderr)
        artifacts_info = None

    return success, artifacts_info, error_msg


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
