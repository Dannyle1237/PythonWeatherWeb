__author__ = 'robswift'

import os
import sys
import itertools
import EB.builder.utilities.csv_interface as csv_interface
import EB.builder.utilities.screener as screener
import EB.builder.utilities.classification as classification
from EB.builder.utilities.ensemble_storage import EnsembleStorage
import EB.builder.utilities.output as output

try:
    import multiprocessing
except ImportError:
    print("multiprocessing requires 2.6 or later")
    sys.exit(1)


def run(itf):
    """
	Run optimize functions.
	"""

    if not itf:
        return 1

    # access user input
    options = SplitInput(itf)

    # read input
    inputpath = os.path.abspath(options.inputpath)
    print(" Reading input file ...")
    molecules = csv_interface.read_csv(inputpath, options)
    if not molecules:
        print("\n '{flag}' was unable to be parsed\n".format(flag=os.path.basename(options.inputpath)))
        sys.exit(1)

    # determine the sort order & ensemble_size
    #sort_order = classification.get_sort_order(molecules)
    sort_order = 'asc'
    ensemble_size = options.ensemble_size

    # loop over all ensembles
    # temp 2/3/15 append to auc_list ef_list & write it out for later histogram construction
    auc_list = []
    ef_list = []
    for size in [x + 1 for x in range(ensemble_size)]:
        auc, ef = optimizor(molecules, sort_order, size, options)
        auc_list += auc
        ef_list += ef
    # temp 2/9/15 write auc_list & ef_list out to files for subsequent post-processing
    f = open('auc_histogram.csv', 'w')
    for value in auc_list:
        f.write('%f\n' % value)
        #f.write('%f, %s\n' % (value[0], value[1]))
    f.close()
    f = open('ef_histogram.csv', 'w')
    for value in ef_list:
        f.write('%f\n' % value)
    f.close()

def optimizor(molecules, sort_order, ensemble_size, options):
    """
	Evaluate the performance of all ensembles of fixed size.
	"""
    # set variables
    ncpu = options.ncpu
    score_field = options.score_field

    # generate an exhaustive list of all possible ensembles
    ensemble_list = make_ensemble_list(molecules, score_field, ensemble_size)

    # set number of processors.
    if not ncpu:
        ncpu = multiprocessing.cpu_count()

    if ncpu > 1:
        print("Determining the performance of {d} ensembles using {n} processors".format(d=len(ensemble_list), n=ncpu))

        if ncpu > len(ensemble_list):
            ncpu = len(ensemble_list)

        jobs = []
        output_queue = multiprocessing.Queue()

        for ensemble_chunk in chunker(ensemble_list, ncpu):
            p = multiprocessing.Process(target=evaluate,
                                        args=(molecules, ensemble_chunk, sort_order, options, output_queue))
            jobs.append(p)
            p.start()

        # collect results into a dictionary
        results = {}
        for i in range(len(jobs)):
            results.update(output_queue.get())

        # stop jobs
        for j in jobs:
            j.join()

    else:
        print("Determining the performance of {d} ensembles using {n} processor".format(d=len(ensemble_list), n=ncpu))
        results = evaluate(molecules, ensemble_list, sort_order, options)

    # peel away the best performing ensemble
    ensemble = screener.find_best_ensemble(results, options)

    # write out the best performing ensemble
    output.write_ensemble(list(ensemble), options)

    # temp 2/9/15 generate and return a list of auc values and ef at fpf = 0.001 to build up a histogram
    nd = max([results[x].ef.keys() for x in results.keys()][0])
    n = int(round(0.001 * nd))
    ef_list = [results[x].get_prop(n, 'ef') for x in results.keys()]
    auc_list = [results[x].get_prop('auc') for x in results.keys()]
    # auc_list = [[results[x].get_prop('auc'), results[x].get_prop('ensemble')] for x in results.keys()]
    return auc_list, ef_list

def evaluate(molecules, ensemble_chunk, sort_order, options, output_queue=None):
    """
	Evaluate VS performance of each ensemble in ensemble_chunk
	"""

    results = {}    # {('receptor_1', ..., 'receptor_n') : ensemble storage object}

    for ensemble in ensemble_chunk:
        results[ensemble] = calculate_performance(molecules, ensemble, sort_order, options)

    if output_queue is not None:
        output_queue.put(results)
    else:
        return results


def calculate_performance(molecules, ensemble, sort_order, options):
    """
    determine the virtual screening performance of the input ensemble, and return the results in an
    ensemble storage object.
    :param molecules:
    :param ensemble:
    :param sort_order: string. either 'asc' (for binding energy estimates)  or 'dsc' (for similarity scores)
    :param options: instance of s
    :return:
    """
    es = EnsembleStorage()
    es.set_prop('ensemble', ensemble)

    # calculate the appropriate score structure type
    score_structure = classification.make_score_structure(molecules, ensemble)

    # determine auc value
    auc_structure = classification.make_auc_structure(score_structure)
    auc = classification.calculate_auc(auc_structure, sort_order, 'no stats')
    es.set_prop('auc', auc)

    # calculate enrichment factors
    for fpf in classification.make_fpfList(options, score_structure):
        fpf = float(fpf)
        ef_structure = classification.make_ef_structure(score_structure, fpf, sort_order)
        if ef_structure:
            ef = classification.calculate_ef(ef_structure, fpf)
            es.set_prop(ef[0], ef[1], 'ef')
    return es


def make_ensemble_list(molecules, score_field, ensemble_size):
    """
	Construct ensemble list
	"""

    # generate list of queries
    queryList = molecules[0].scores.keys()

    # nchoosek
    ensemble_iterator = itertools.combinations(queryList, ensemble_size)

    # list of tuples: [(query1, query2), ... (queryN-1, queryN)
    ensembleList = []

    for ensemble in ensemble_iterator:
        ensembleList.append(ensemble)

    return ensembleList


def chunker(ensemble_list, ncpu):
    """
	Generate successive chunks of ensemble_list.
	"""

    # determine sublist lengths
    length = int(len(ensemble_list) / ncpu)

    # generator
    for i in range(0, len(ensemble_list), length):
        yield ensemble_list[i:i + length]


class SplitInput:
    """
	Return command line input or print error message and return 1.
	"""

    def __init__(self, itf):
        self.inputpath = itf.get_string('--input')
        self.framework = itf.get_string('--framework')
        self.outname = itf.get_string('--outname')
        self.ensemble_size = itf.get_int('--ensemble_size')
        self.fpf = itf.get_float('--fpf')
        self.ncpu = itf.get_int('--ncpu')
        self.score_field = 'receptor'
        self.status_field = 'status'
        self.active_label = '1'
        self.decoy_label = '0'