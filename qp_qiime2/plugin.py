# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

import os

from biom import example_table
from qiime.sdk import PluginManager, Artifact, SubprocessExecutor

def job_id_to_plugin_workflow_name(job_id):
    # the available QIIME 2 plugins and their associated workflows should
    # be made available through another entry point, so Qiita can ask one of
    # its plugins what workflows it has available
    # TODO: actually pull this info from Qiita
    return ('feature-table', 'rarefy')

def job_id_to_input_artifacts(job_id, artifact_names):
    # TODO: actually pull this info from Qiita
    return [('table', example_table)]

def job_id_to_input_parameters(job_id, parameter_names):
    # TODO: actually pull this info from Qiita
    return {'depth': 2}


def job_id_to_workflow_executor_arguments(job_id, plugin_manager):
    plugin_name, workflow_name = job_id_to_plugin_workflow_name(job_id)
    wf = plugin_manager.plugins[plugin_name].workflows[workflow_name]

    input_artifact_fps = {}
    for artifact_name, artifact_data in job_id_to_input_artifacts(job_id):
        # TODO: actually make this a randomly generated file name
        temp_file_name = "random-in.qtf"
        input_artifact_fps[artifact_name] = temp_file_name
        Artifact.save(artifact_data,
                      wf.signature.input_artifacts[artifact_name],
                      None,
                      temp_file_name)
    input_parameters = job_id_to_input_parameters(job_id)
    # TODO: actually make this a randomly generated file name
    output_artifacts = {name: 'random-out.qtf'
                        for name in wf.signature.output_artifacts.keys()}
    return wf, input_artifact_fps, input_parameters, output_artifacts


def execute_job(url, job_id, output_dir):
    plugin_manager = PluginManager()

    wf, input_artifact_fps, input_parameters, output_artifacts = \
        job_id_to_workflow_executor_arguments(job_id, plugin_manager)

    # execute workflow
    executor = SubprocessExecutor()
    future_ = executor(wf,
                       input_artifact_fps,
                       input_parameters,
                       output_artifacts)

    completed_process = future_.result()
    if completed_process.returncode == 0:
        # send data to Qiita
        print(output_artifacts)
    else:
        print(completed_process.stdout)
        exit(-1)
