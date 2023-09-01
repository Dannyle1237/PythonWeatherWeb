__author__ = 'robswift'

import csv
import os
import re

def write_roc(roc_structure, inputfilename, options, fw_type = None):
    """
    writes ROC output
    :param roc_structure: a.k.a auc_structure, generated in /common_tools/classification/make_auc_structure()
    :param inputfilename: name of the file specified on the command line that contains the ensemble. This will be
    used to generate the name of the file that contains the ROC data.
    :param options:
    :return:
    """
    # check if ./ROC_DATA exists. If not, create it
    rocdir = os.path.join(os.getcwd(), 'ROC_DATA')
    if not os.path.exists(rocdir):
        os.makedirs(rocdir)

    # generate the root name of the file that will hold the roc data
    if len(inputfilename.split('.')) == 1:
        rootname = inputfilename.split('.')
    else:
        rootname = inputfilename.split('.csv')[0]

    # add '.csv' to the root name to give the name of the file
    if fw_type:
        filename = '%s_fw_roc.csv' % rootname
    else:
        filename = '%s_roc.csv' % rootname

    # the path to the file
    file = os.path.join(rocdir, filename)

    # open file and create a csv writer object
    f = open(file, 'w')
    rocwriter = csv.writer(f)

    # create header
    header = ['id', 'score', 'score source', 'status', 'fpf', 'tpf']
    rocwriter.writerow(header)

    # write contents
    for tup in roc_structure:
        rocwriter.writerow(list(tup))

    f.close()

def write_diff_summary(stats, options):
    """
    writes contents of stats to outname_diff_stats.csv
    :param stats:
    :param options:
    """
    # set output file name
    out_name = '%s_diff_stats.csv' % options.outname
    filename = os.path.join(os.getcwd(), out_name)

    # open file & create csv writer object
    f = open(filename, 'w')
    statswriter = csv.writer(f)

    # write header
    statswriter.writerow(stats['header'])

    # write auc
    line = ['AUC']
    for metric_tup in stats['AUC']:
        try:
            len(metric_tup)
            line.append('[%0.2f, %0.2f]' % (metric_tup[0], metric_tup[1]))
        except:
            line.append('%0.6f' % metric_tup)
    statswriter.writerow(line)

    # write enrichment factors
    for title in [x for x in sorted(stats.keys()) if 'E' in x]:
        line = [title]
        for metric_tup in stats[title]:
            try:
                len(metric_tup)
                line.append('[%0.2f, %0.2f]' % (metric_tup[0], metric_tup[1]))
            except:
                line.append('%0.12f' % metric_tup)
        statswriter.writerow(line)

def write_summary(stats, options, fw_type = None):
    """
    writes contents of stats to outname_stats.csv
    :param stats: {'num_of_queries' : [(auc, auclow, auchigh),(fpf, ef, eflow, efhigh),(fpf, ef, eflow, efhigh),...,]}
    if fw_type is not None, stats has the form given below
    {'num_of_queries' : [(auc, auclow, auchigh), (fpf, ef, eflow, efhigh), (fpf, acd), ...,], }
    :param options: postanalyze.SplitInput object
    """

    # set output file name
    out_name = '%s_stats.csv' % options.outname
    filename = os.path.join(os.getcwd(), out_name)

    # open file & create csv writer object
    f = open(filename, 'w')
    stats_writer = csv.writer(f)

    # sort the stats keys using the "decorate, sort, un decorate method"
    decorated = [(x, int(re.search(r'_([0-9]|[0-9][0-9]|[0-9][0-9][0-9])_queries.csv$', x, re.M).group().split('_')[1]))
                 for x in stats.keys() if '_queries.csv' in x]
    decorated.sort(key=lambda element: element[1])
    sorted_stats_keys = [x[0] for x in decorated]

    for leftovers in [x for x in stats.keys() if '_queries.csv' not in x]:
        sorted_stats_keys.append(leftovers)

    for ensemble in sorted_stats_keys:
        # write header & first line
        if f.tell() == 0:
            header = ['Ensemble name', 'AUC', '95% CI']
            line = ['%s' % ensemble]
            for metric_tup in stats[ensemble]:
                if len(metric_tup) == 3:
                    line.append('%0.4f' % metric_tup[0])
                    line.append('[%0.2f, %0.2f]' % (metric_tup[1], metric_tup[2]))
                if len(metric_tup) == 4:
                    header.append('E%s' % metric_tup[0])
                    header.append('95% CI')
                    line.append('%0.2f' % metric_tup[1])
                    line.append('[%0.2f, %0.2f]' % (metric_tup[2], metric_tup[3]))
            stats_writer.writerow(header)
            stats_writer.writerow(line)
        # write remaining lines
        else:
            line = ['%s' % ensemble]
            for metric_tup in stats[ensemble]:
                if len(metric_tup) == 3:
                    line.append('%0.4f' % metric_tup[0])
                    line.append('[%0.2f, %0.2f]' % (metric_tup[1], metric_tup[2]))
                if len(metric_tup) == 4:
                    line.append('%0.4f' % metric_tup[1])
                    line.append('[%0.2f, %0.2f]' % (metric_tup[2], metric_tup[3]))
            stats_writer.writerow(line)

    # write maximum values
    line = ['max values:']
    for metric_tup in stats[list(stats.keys())[0]]:
        if len(metric_tup) == 3:
            line.append('1.00') # auc
            line.append(' ') # CI space
        if len(metric_tup) == 4:
            line.append('%0.2f' % (1 / metric_tup[0])) # max ef
            line.append(' ') # CI space
    stats_writer.writerow(line)

    # close file
    f.close()


def write_ensemble(ensemble, options):
    """
	Prints out the ensemble composition at each size
	"""

    # set output file name
    size = len(ensemble)
    filename = '%s_%s_queries.csv' % (options.outname, size)
    file = os.path.join(os.getcwd(), filename)

    f = open(file, 'w')

    out = ', '.join(ensemble)

    f.write(out)

    f.close()