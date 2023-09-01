__author__ = 'robswift'

import os
import sys
from EB.builder.utilities import classification
from EB.builder.utilities.ensemble_storage import  EnsembleStorage
from EB.builder.utilities import csv_interface
from EB.builder.utilities import output
from EB.builder.utilities import screener

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
        print("\n '{f} was unable to be parsed\n".format(f=os.path.basename(options.inputpath)))
        sys.exit(1)
    approximator(molecules, options)

def approximator(molecules, options, sort_order=None):
    ensemble_size = options.ensemble_size

    if not sort_order:
        sort_order = classification.get_sort_order(molecules)

    print("Performing calculations")
    results = rank_queries(molecules, sort_order, options)

    ensemble = []
    while results:
        ensemble = construct_ensemble(results, ensemble, options)
        output.write_ensemble(ensemble, options)
        if len(ensemble) == ensemble_size:
            return 0

def construct_ensemble(results, ensemble, options):
    """

    :param results: {'(query)' : EB.builder.utilities.ensemble_storage.EnsembleStorage object}
    :param ensemble:
    :param options:
    :return:
    """
    best_query = screener.find_best_ensemble(results, options)
    ensemble.append(best_query[0])
    del results[best_query]
    return list(ensemble)

def rank_queries(molecules, sort_order, options):
    results = {}
    for query in [query for query in list(molecules[0].scores.keys())]:
        formatted_query = []
        formatted_query.append(query)
        formatted_query = tuple(formatted_query)
        es = EnsembleStorage()
        es.set_prop('ensemble', formatted_query)
        score_structure = classification.make_score_structure(molecules, formatted_query)
        auc_structure = classification.make_auc_structure(score_structure)
        auc = classification.calculate_auc(auc_structure, sort_order, 'no stats')
        es.set_prop('auc', auc)
        for fpf in classification.make_fpfList(options, score_structure):
            fpf = float(fpf)
            ef_structure = classification.make_ef_structure(score_structure, fpf, sort_order)
            if ef_structure:
                ef = classification.calculate_ef(ef_structure, fpf)
                es.set_prop(ef[0], ef[1], 'ef')
        results[formatted_query] = es
    return results

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
