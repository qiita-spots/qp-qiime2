-- Add the plugin
INSERT INTO qiita.software (name, version, description, environment_script, start_script, software_type_id)
    VALUES ('QIIME2', '2.0.0', 'Qiime 2', 'workon qp-qiime2', 'start-qiime2', 1);

-- Add the commands
INSERT INTO qiita.software_command (software_id, name, description) VALUES
    (4, 'feature-table: rarefy', 'Rarefies an OTU table');

-- Add the parameters
INSERT INTO qiita.command_parameter (command_id, parameter_name, parameter_type, required)
    VALUES (8, 'table', 'artifact', True),
           (8, 'depth', 'integer', False);

-- Link the parameter type with the artifact type
INSERT INTO qiita.parameter_artifact_type (command_parameter_id, artifact_type_id)
    VALUES (43, 7);

-- Add the output
INSERT INTO qiita.command_output (name, command_id, artifact_type_id)
    VALUES ('rarefied_table', 8, 7);

-- Add a default parameter set
INSERT INTO qiita.default_parameter_set (command_id, parameter_set_name, parameter_set)
    VALUES (8, 'test-rarefaction', '{"depth":100}');

INSERT INTO qiita.oauth_software (client_id, software_id)
    VALUES ('4MOBzUBHBtUmwhaC258H7PS0rBBLyGQrVxGPgc9g305bvVhf6h', 4);
