__author__ = 'robswift'


import os
import sys
import EB.builder.utilities.csv_interface as csv_interface
import EB.builder.utilities.screener as screener
import EB.builder.utilities.classification as classification
from EB.builder.utilities.ensemble_storage import EnsembleStorage
import EB.builder.utilities.output as output


def run(itf):
    """
    run approximate functions
    :param itf:
    :return:
    """

    if not itf:
        return 1

    # access user input
    options = SplitInput(itf)

    # read input
    print(" Reading input file ...")
    molecules = csv_interface.read_csv(os.path.abspath(options.inputpath), options)
    if not molecules:
        print("\n '{f}' was unable to be parsed\n".format(f=os.path.abspath(options.inputpath)))
        sys.exit(1)

    if approximator(molecules, options) == 1:
        return 0


def approximator(molecules, options, sort_order=None, frameworks=[], ensemble=[]):
    """
    recursively rank queries
    :param molecules:
    :param options:
    :param sort_order:
    :param frameworks:
    :param ensemble:
    :return:
    """

    # set variables
    ensemble_size = options.ensemble_size

    if not sort_order:
        sort_order = classification.get_sort_order(molecules)

    # construct ensemble
    print("Performing calculations for ensemble size {s}".format(s=(len(ensemble) + 1)))
    ensemble = rank_queries(molecules, ensemble, sort_order, options)

    # write stats & ensemble
    output.write_ensemble(list(ensemble), options)

    if len(ensemble) == ensemble_size:
        return 1
    else:
        return approximator(molecules, options, sort_order, frameworks, ensemble)


def rank_queries(molecules, ensemble, sort_order, options):
    """
    rank queries by value added to existing ensemble
    :param molecules:
    :param score_field:
    :param ensemble:
    :param sort_order:
    :param options:
    :return:
    """

    # generate query list
    query_list = [x for x in list(molecules[0].scores.keys()) if x not in ensemble]

    results = {}

    for query in query_list:
        es = EnsembleStorage() # an ensemble storage project

        # generate test_ensemble
        test_ensemble = ensemble[0:]
        test_ensemble.append(query)
        test_ensemble = tuple(test_ensemble)
        es.set_prop('ensemble', test_ensemble)

        # calculate its performance
        score_structure = classification.make_score_structure(molecules, test_ensemble)

        # determine auc value
        auc_structure = classification.make_auc_structure(score_structure)
        auc = classification.calculate_auc(auc_structure, sort_order, 'no stats')
        es.set_prop('auc', auc)

        # if the enrichment factor was set to anything other than 1, then we're training to maximize the corresponding
        # enrichment factor
        for fpf in classification.make_fpfList(options, score_structure):
            fpf = float(fpf)
            ef_structure = classification.make_ef_structure(score_structure, fpf, sort_order)
            if ef_structure:
                ef = classification.calculate_ef(ef_structure, fpf)
                es.set_prop(ef[0], ef[1], 'ef')

        # append results to metric list
        results[test_ensemble] = es

    # peel away the best performing ensemble
    best_ensemble = screener.find_best_ensemble(results, options)

    return list(best_ensemble)


class SplitInput:
    """
    interfaces with an instance of the input interface class to return command line arguments
    :param itf: input interface object
    """

    def __init__(self, itf):
        self.inputpath = itf.get_string('--input')
        self.framework = itf.get_string('--framework')
        self.outname = itf.get_string('--outname')
        self.ensemble_size = itf.get_int('--ensemble_size')
        self.fpf = itf.get_float('--fpf')
        self.score_field = 'receptor'
        self.status_field = 'status'
        self.active_label = '1'
        self.decoy_label = '0'