__author__ = 'robswift'

import os
import gzip
import sys
import csv
import re
from EB.builder.utilities.molecules import Molecule


def read_csv(csv_file, options, ensemble_list=None):
    """
	Read csv and return molList, otherwise print error and exit.
	"""
    name, ext = os.path.splitext(csv_file)
    try:
        if ext == '.gz':
            f = gzip.open(csv_file, 'rb')
        else:
            f = open(csv_file, 'rU')
    except IOError:
        print(" \n '{f}' could not be opened\n".format(f=os.path.basename(csv_file)))
        sys.exit(1)


    csv_reader = csv.reader(f)
    molList = []
    line_number = 1

    for line in csv_reader:

        if line_number == 1:
            if ensemble_list:
                prop_indices = read_header(line, options, ensemble_list)
            else:
                prop_indices = read_header(line, options)
        else:
            mol = Molecule()
            if ensemble_list:
                mol = read_line(line, options, prop_indices, mol, ensemble_list)
            else:
                mol = read_line(line, options, prop_indices, mol)
            if mol == 1:
                print(" skipping molecule {m}\n".format(m=(line_number - 1)))
            else:
                molList.append(mol)
        line_number += 1

    return molList


def read_header(header_labels, options, ensemble_list=None):
    """
	read csv header, returns prop_indices dictionary, {'label' : index}
	where label is the csv column label, and index the column index
	"""

    # set input variables
    filename = os.path.basename(options.inputpath)
    if ensemble_list:
        score_field = ensemble_list[0]
        ensemble_size = len(ensemble_list)
    else:
        score_field = options.score_field
        ensemble_size = options.ensemble_size


    status_field = options.status_field

    # check on score field values
    if ensemble_list:
        query_list = ensemble_list
        # verify that the input queries are in the inputfille
        if not [x for x in query_list if x in header_labels]:
            out = ' '.join(query_list)
            print("\n '{q}' are not in '{e}'\n".format(q=out, e=os.path.basename(filename)))
            sys.exit(1)
    else:
        score_field_matcher = re.compile(score_field)
        query_list = [x for x in header_labels if score_field_matcher.match(x)]

        if not query_list:
            print("\n There are no '{s}' columns in '{c}'\n".format(s=score_field, c=os.path.basename(filename)))
            sys.exit(1)

        # confirm that there are a sufficient no. of queries to generate
        # the requested ensemble size
        if len(query_list) < ensemble_size:
            print("\n '{s}' will support a maximum ensemble size of '{d}'\n".format(s=os.path.basename(filename),
                                                                              d=len(query_list)))
            sys.exit(1)

    # confirm that status field exists
    status_field_matcher = re.compile(status_field)
    if not any(x for x in header_labels if status_field_matcher.match(x)):
        print("\n The status_field column '{s}' doesn't exist\n".format(s=status_field))
        sys.exit(1)

    # initiate prop_column dictionary and csv column index
    prop_indices = {}
    index = 0

    for label in header_labels:
        prop_indices[label] = index
        index += 1

    return prop_indices


def read_line(csv_contents, options, prop_indices, mol, ensemble_list=None):
    """
	read csv line
	"""
    if not ensemble_list:
        score_field = options.score_field
    status_field = options.status_field
    active_label = options.active_label
    decoy_label = options.decoy_label

    # do the active/decoy labels have appropriate values?
    active_value_matcher = re.compile(active_label)
    decoy_value_matcher = re.compile(decoy_label)
    status_label_index = prop_indices[status_field]

    if not active_value_matcher.match(csv_contents[status_label_index]) and not decoy_value_matcher.match(
            csv_contents[status_label_index]):
        print("\n molecule lacks appropriate status label")
        return 1

    # are the score field values defined?
    score_field_indices = []
    if ensemble_list:
        queryList = ensemble_list
    else:
        queryList = [x for x in prop_indices.keys() if score_field in x]
    for query in queryList:
        score_field_indices.append(prop_indices[query])
    for value in [csv_contents[x] for x in score_field_indices]:
        if value in ('', 'n/a', 'N/A', None):
            print("\n molecule lacks appropriate score field value")
            return 1

    # loop over property values
    for label in prop_indices.keys():

        # get property value
        value_index = prop_indices[label]
        value = csv_contents[value_index]

        # set corresponding molecule attribute
        if label in queryList:
            mol.SetProp(label, value, 'score')
        else:
            mol.SetProp(label, value)

    # return mol
    return mol