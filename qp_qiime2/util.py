# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------


def get_qiime2_type_name_and_predicate(element):
    """helper method to get the qiime2 type name and predicate

    Parameters
    ----------
    element : qiime2.core.type.signature
        The signature to parse

    Returns
    -------
    str, dict
        The name and the predicate of the inputed signature
    """
    to_ast = element.qiime_type.to_ast()
    if to_ast['type'] == 'union':
        # union types allow to give choices to another type of parameter; for
        # example for a range(1, 10) give the choice `ignore`. These are not
        # necessary in Qiita as they are simply ignored if unchanged. Thus,
        # we loop over the members of the union and ingore `Choices`.
        to_ast = [x for x in to_ast['members']
                  if x['predicate'] is None or
                  x['predicate']['name'] != 'Choices'][0]
        predicate = to_ast['predicate']
    elif to_ast['name'] == 'FeatureData':
        predicate = []
        for f in to_ast['fields']:
            if 'members' in f:
                for fm in f['members']:
                    predicate.append(fm['name'])
            elif 'mapping' in f:
                for fm in f['mapping']:
                    for fme in fm:
                        predicate.append(fme['name'])
            else:
                predicate.append(f['name'])
        predicate = sorted(list(set(predicate)))
    elif to_ast['name'] == 'List' and element.qiime_type.predicate is None:
        # just taking the first one to keep things rolling, required by
        # PluginManager().plugins['composition'].methods[
        #     'ancombc'].signature.parameters['reference_levels']
        predicate = to_ast['fields'][0]['name']
    else:
        predicate = element.qiime_type.predicate
    name = to_ast['name']

    return name, predicate
