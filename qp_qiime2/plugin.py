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
from os import environ
from configparser import ConfigParser

from qiita_client import QiitaClient

from biom import example_table
from qiime.sdk import PluginManager, Artifact, SubprocessExecutor


def get_plugin_workflow_lookup():
    plugin_manager = PluginManager()
    result = {}
    for plugin_name, plugin in plugin_manager.plugins.items():
        for workflow_name, workflow in plugin.workflows:
            result['%s: %s' % (plugin_name, workflow_name)] = workflow
    return result


def retrieve_artifact_info(artifact_id, artifact_type):
    # TODO: retrieve the artifact information from the Qiita server
    artifact_data = example_table
    # TODO: actually make this a randomly generated file name
    temp_file_name = "random-in.qtf"
    # TODO: None is the provenance
    Artifact.save(artifact_data, artifact_type, None, temp_file_name)
    return temp_file_name


def job_id_to_workflow_executor_arguments(job_id, job_info):
    wf_lookup = get_plugin_workflow_lookup()
    wf = wf_lookup[job_info['command']]
    parameters = job_info['parameters']

    input_artifact_fps = {}
    for ia_name, ia_type in wf.signature.input_artifacts.items():
        artifact_fp = retrieve_artifact_info(parameters[ia_name])
        input_artifact_fps[ia_name] = artifact_fp

    input_parameters = {ip_name: parameters[ip_name]
                        for ip_name in wf.signature.input_parameters}

    # TODO: actually make this a randomly generated file name
    output_artifacts = {name: 'random-out.qtf'
                        for name in wf.signature.output_artifacts.keys()}

    return wf, input_artifact_fps, input_parameters, output_artifacts


def format_artifacts_info(output_artifacts):
    pass


def execute_task(job_id, job_info):
    wf, input_artifact_fps, input_parameters, output_artifacts = \
        job_id_to_workflow_executor_arguments(job_id, job_info)

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
        success, artifacts_info, error_msg = execute_task(job_id, job_info)
    except Exception:
        exc_str = repr(traceback.format_exception(*sys.exc_info()))
        error_msg = ("Error executing %s:\n%s" % (job_info['command'],
                                                  exc_str))
        success = False
        artifacts_info = None

    # The job completed
    qclient.complete_job(job_id, success, error_msg=error_msg,
                         artifacts_info=artifacts_info)
