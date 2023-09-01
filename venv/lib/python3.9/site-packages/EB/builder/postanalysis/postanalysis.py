__author__ = 'robswift'

import os
import sys
import copy
import EB.builder.utilities.csv_interface as csv_interface
import EB.builder.utilities.classification as classification
import EB.builder.utilities.output as output


def run(itf):
    """
	Run postanalyze functions.
	"""

    if not itf:
        return 1

    # access user input
    options = SplitInput(itf)

    # check input args
    error_check(options)

    # read input files
    try:
        molecules, ensemble_lookup = ReadFiles(options)
    except:
        return 1

    if options.compare:
        compare(molecules, ensemble_lookup, options)
    else:
        evaluate_list(molecules, ensemble_lookup, options)


def compare(molecules, ensemble_lookup, options):
    """
    compare stuff
    :param molecules:
    :param ensemble_lookup:
    :param options:
    :return:
    """

    print(" Analyzing differences ... ")
    print('')
    sort_order = classification.get_sort_order(molecules)

    ensemble1 = sorted(ensemble_lookup.keys())[0]
    ensemble2 = sorted(ensemble_lookup.keys())[1]

    stats = {}
    stats['header'] = [' ']
    name = os.path.basename(ensemble1).replace('.csv', '')
    stats['header'].append(name)
    name = os.path.basename(ensemble2).replace('.csv', '')
    stats['header'].append(name)
    stats['header'].append('Difference')
    stats['header'].append('95% CI')
    stats['header'].append('p-value')

    molecules1 = copy.deepcopy(molecules)
    molecules2 = copy.deepcopy(molecules)

    score_structure1 = classification.make_score_structure(molecules1, ensemble_lookup[ensemble1])
    score_structure2 = classification.make_score_structure(molecules2, ensemble_lookup[ensemble2])

    auc_structure_1 = classification.make_auc_structure(score_structure1)
    auc_structure_2 = classification.make_auc_structure(score_structure2)

    # calculate auc value differences
    auc_diff = classification.calculate_auc_diff(auc_structure_1, auc_structure_2, sort_order)

    stats['AUC'] = auc_diff

    # calculate enrichment factor differences
    fpfList = make_fpfList(options)
    for fpf in fpfList:
        fpf = float(fpf)
        ef_structure1 = classification.make_ef_structure(score_structure1, fpf, sort_order)
        ef_structure2 = classification.make_ef_structure(score_structure2, fpf, sort_order)

        if ef_structure1 and ef_structure2:
            ef_diff = classification.calculate_ef_diff(ef_structure1, ef_structure2, fpf)
            title = 'E%s' % fpf
            stats[title] = ef_diff

    # write results summary
    output.write_diff_summary(stats, options)

    # write roc curves
    if options.write_roc:
        print(" Writing ROC data ... ")
        print('')
        output.write_roc(auc_structure_1, ensemble1, options)
        output.write_roc(auc_structure_2, ensemble2, options)

    # plot
    if options.plot:
        print(" Making plots ... ")
        print('')
        plotter(molecules, ensemble_lookup, options)


def evaluate_list(molecules, ensemble_lookup, options):
    """
    Evaluate a list of ensembles and return statistics and ROC plots if appropriate
    """

    # create stats dictionaries to store results from each ensemble
    stats = {}  # {file name : metric_List}

    # print progress messages
    if options.write_roc:
        print(" Determining virtual screening performance and writing ROC data ... ")
        print('')
    else:
        print(" Determining virtual screening performance ...")
        print('')

    for filename in sorted(ensemble_lookup.keys()):
        metric_List = calculate_metrics(molecules, ensemble_lookup, filename, options)
        stats[filename] = metric_List

    # write results summary
    output.write_summary(stats, options, fw_type = None)

    # plot
    if options.plot:
        print(" Making plots ... ")
        print
        plotter(molecules, ensemble_lookup, options)


def calculate_metrics(molecules, ensemble_lookup, filename, options):
    """
    Determine the virtual screening performance of the ensemble
    :param molecules: list [mol_object_1, mol_object_2, .... ]
    :param ensemble: tuple (receptor_x, receptor_y, .... )
    :param options: interface object that makes command line arguments available.
    :return:
    """
    metric_List = []    # [(auc, auclow, auchigh), (fpf, ef, eflow, efhigh), (fpf, ef, eflow, efhigh), ..., ]
    sort_order = 'asc'

    # set up the appropriate score_structure data
    score_structure = classification.make_score_structure(molecules, ensemble_lookup[filename])

    # calculate auc values
    auc_structure = classification.make_auc_structure(score_structure)
    auc = classification.calculate_auc(auc_structure, sort_order)
    metric_List.append(auc)

    # calculate enrichment factor values
    for fpf in make_fpfList(options):
        fpf = float(fpf)
        ef_structure = classification.make_ef_structure(score_structure, fpf, sort_order)
        if ef_structure:
            ef = classification.calculate_ef(ef_structure, fpf, None, 'include_intervals')
            metric_List.append(ef)

    if options.write_roc:
        output.write_roc(auc_structure, filename, options)

    return metric_List


def make_fpfList(options):
    """
    aggregate default fpf values and user-defined fpf values where enrichment factor calculations will be attempted.
    :param options: SplitInput object
    :return defaults: list includes default fpf values & unique user-defined values
    """

    user_values = options.fpf

    defaults = ['0.0001', '0.001', '0.01', '0.05']

    if user_values:
        for fpf in user_values:
            if fpf not in defaults:
                defaults.append(fpf)

        defaults.sort()

    return defaults


def error_check(options):
    """
    Error check
    :rtype : object
    """

    compare = options.compare
    ensemble_paths = options.ensemble_paths

    if compare and len(ensemble_paths) > 2:
        print("\n Only 2 ensembles can be compared, {d} were specified\n".format(d=len(ensemble_paths)))
        sys.exit(1)


def plotter(molecules, ensemble_lookup, options):
    """
    plot ROC curves for ensembles in ensemble_lookup
    :param molecules:
    :param ensemble_lookup:
    :param options:
    :return:
    """

    try:
        import matplotlib

        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n Plotting requires matplotlib to be installed\n")
        sys.exit(1)

    for ensemble in ensemble_lookup.keys():

        # create figure
        fig = plt.figure()

        # create the queries subplot, the left subplot
        # create the left hand subplot
        ax1 = fig.add_subplot(121)

        for query in sorted(ensemble_lookup[ensemble]):
            query_list = []
            query_list.append(query)
            score_structure = classification.make_score_structure(molecules, query_list)
            auc_structure = classification.make_auc_structure(score_structure)
            tpf = []
            fpf = []
            for mol in auc_structure:
                fpf.append(mol[4])
                tpf.append(mol[5])

            # add axis-labels and a title
            ax1.set_xlabel('FPF')
            ax1.set_ylabel('TPF')
            title = 'query performance'
            ax1.set_title(title)

            # add plot data and labels for the legend
            lbl = query
            ax1.plot(fpf, tpf, lw=3, label=lbl)

        # get legend handles and labels, then reverse their order
        handles, labels = ax1.get_legend_handles_labels()
        ax1.legend(handles[::-1], labels[::-1])

        # add the legend
        ax1.legend(handles, labels, loc='best')

        # create the ensemble subplot, the right subplot
        score_structure = classification.make_score_structure(molecules, ensemble_lookup[ensemble])
        auc_structure = classification.make_auc_structure(score_structure)
        tpf = []
        fpf = []
        for mol in auc_structure:
            fpf.append(mol[4])
            tpf.append(mol[5])

        # create right hand subplot
        ax2 = fig.add_subplot(122)

        # add axis-labels and a title
        ax2.set_xlabel('FPF')
        ax2.set_ylabel('TPF')
        title = 'ensemble performance'
        ax2.set_title(title)

        # add plot data and a label for the legend
        lbl = 'ensemble'
        ax2.plot(fpf, tpf, lw=3, label=lbl)

        # get legend handles and labels, then reverse their order
        handles, labels = ax2.get_legend_handles_labels()
        ax2.legend(handles[::-1], labels[::-1])

        # add the legend
        ax2.legend(handles, labels, loc='best')

        # save figure
        figurename = options.outname + '_' + ensemble.replace('.csv', '') + '.pdf'
        filename = os.path.join(os.getcwd(), figurename)
        plt.savefig(filename, bbobx='tight', format='pdf')


def ReadFiles(options):
    # set variables
    ensemble_paths = options.ensemble_paths

    # read ensembles into ensemble lookup {'Ensemble_N' : [query_name_list]}
    ensemble_lookup = {}

    for ensemble_path in ensemble_paths:
        try:
            ensemble_file = open(ensemble_path, 'r')
        except IOError:
            print("\nUnable to open ensemble_list: {l}\n".format(l=ensemble_path))
            return 1

        ensemble_queries = [query.strip() for query in ensemble_file.read().split(',')]
        ensemble_file.close()
        if len(ensemble_queries) == 0 or '' in ensemble_queries:
            print("\n{l} is empty\n".format(l=ensemble_path))
            return 1

        key = os.path.basename(ensemble_path)
        ensemble_lookup[key] = ensemble_queries

    # Run consistency checks on the input csv
    uniq = []
    for ensemble in ensemble_lookup.keys():
        for unique_query in [query for query in ensemble_lookup[ensemble] if query not in uniq]:
            uniq.append(unique_query)

    # read input csv
    inputpath = os.path.abspath(options.inputpath)
    print('')
    print(" Reading input file ...")
    print('')
    molecules = csv_interface.read_csv(inputpath, options, uniq)
    if not molecules:
        print("\n '%s' was unable to be parsed\n")
        sys.exit(1)

    return molecules, ensemble_lookup


class SplitInput:
    """
	Return command line input or print error message and return 1.
	"""

    def __init__(self, itf):
        self.inputpath = itf.get_string('--input')
        self.outname = itf.get_string('--outname')
        self.fpf = itf.get_list('--fpf')
        self.plot = itf.get_boolean('--plot')
        self.write_roc = itf.get_boolean('--write_roc')
        self.compare = itf.get_boolean('--compare')
        self.ensemble_paths = itf.get_list('--ensemble_list')
        self.status_field = 'status'
        self.active_label = '1'
        self.decoy_label = '0'