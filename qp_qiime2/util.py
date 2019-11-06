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
        # union types allow to give choices to another type of paramter; for
        # example for a range(1, 10) give the choice `ignore`. These are not
        # necessary in Qiita as they are simply ignored if unchanged. Thus,
        # we loop over the members of the union and ingore `Choices`.
        to_ast = [x for x in to_ast['members']
                  if x['predicate']['name'] != 'Choices'][0]
        predicate = to_ast['predicate']
    else:
        predicate = element.qiime_type.predicate
    name = to_ast['name']

    return name, predicate
