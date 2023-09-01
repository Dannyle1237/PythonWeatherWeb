__author__ = 'robswift'

import os
import sys
import re
import csv
import gzip
import random


def run(itf):
    """
	Run preprocess functions
	"""

    if not itf:
        return 1

    # access command-line arguments
    options = SplitInput(itf)

    # read input
    infile = os.path.abspath(options.input)
    molList = read_csv(infile, options)

    # split molList into actives and decoys
    activeList, decoyList = partition(molList, options)

    # split actives and decoys into training and validation sets
    trainset, valset = split(activeList, decoyList, options)

    # write csv files formatted for ensemble builder
    csv_writer(trainset, options, 'training_set')
    csv_writer(valset, options, 'test_set')


def write_csv_header(mol, csv_writer):
    """
	Write the csv header
	"""

    # create line list where line elements for writing will be stored
    line = []

    # ID
    line.append('id')

    # status
    line.append('status')

    # query labels
    queryList = mol.properties.keys()
    for queryLabel in queryList:
        line.append(queryLabel)

    # write line
    csv_writer.writerow(line)


def write_csv_line(mol, csv_writer, options):
    """
	Parse mol object and write a line to the csv file
	"""

    # set variables
    status_field = options.status_field

    # elements for writing will be stored in the line list
    line = []

    # ID
    id = mol.GetProp('id')
    if id is not None:
        line.append(id)
    else:
        line.append('n/a')

    # status
    line.append(mol.GetProp(status_field))

    # query labels
    queryList = mol.properties.keys()
    for queryLabel in queryList:
        line.append(mol.properties[queryLabel])

    # write line
    csv_writer.writerow(line)


def csv_writer(molecules, options, prefix):
    """
	Write a csv file.
	"""

    # output file
    outdir = os.getcwd()
    filename = prefix + '.csv'
    outfile = os.path.join(outdir, filename)

    # initiate csv writer object
    f = open(outfile, 'w')
    csv_writer = csv.writer(f)

    # write csv header
    mol = molecules[0]
    write_csv_header(mol, csv_writer)

    # write csv lines
    for mol in molecules:
        write_csv_line(mol, csv_writer, options)

    # close file
    f.close()


def split(activeList, decoyList, options):
    """
	Create training and validation sets
	"""

    # set input variables
    training_fraction = options.training_fraction
    decoy_to_active = options.decoy_to_active

    # take care of default decoy_to_active ratio
    if decoy_to_active is None:
        decoy_to_active = len(decoyList) / len(activeList)

    # verify that there are enough molecules to satisfy the ratio
    if len(decoyList) < (len(activeList) * decoy_to_active):
        max = len(decoyList) / len(activeList)
        print("\n The maximum decoy to active ratio the input file will support is {f} \n".format(f=max))
        sys.exit(1)

    # randomly split the actives
    trainsize = int(round(training_fraction * len(activeList)))

    trainIndex = []
    valIndex = []
    trainIndex = random.sample(range(len(activeList)), trainsize)
    valIndex = [x for x in range(len(activeList)) if x not in trainIndex]

    trainactives = [activeList[index] for index in trainIndex]
    valactives = [activeList[index] for index in valIndex]

    # match up decoys
    trainsize = len(trainactives) * decoy_to_active
    valsize = len(valactives) * decoy_to_active

    trainIndex = []
    valIndex = []
    trainIndex = random.sample(range(len(decoyList)), int(trainsize))
    valIndex = [x for x in range(len(decoyList)) if x not in trainIndex][0:int(valsize)]

    traindecoys = [decoyList[index] for index in trainIndex]
    valdecoys = [decoyList[index] for index in valIndex]

    # merge actives and decoys for each set
    trainset = trainactives + traindecoys
    valset = valactives + valdecoys

    # return sets
    return trainset, valset


def partition(molList, options):
    """
	Partition molList into activeList and decoyList
	"""
    # set input variables
    status_field = options.status_field
    active_label = options.active_label
    decoy_label = options.decoy_label

    # initiate lists
    activeList = []
    decoyList = []

    # partition moList
    for mol in molList:
        if mol.GetProp(status_field) == active_label:
            activeList.append(mol)
        elif mol.GetProp(status_field) == decoy_label:
            decoyList.append(mol)

    # return partitions
    return activeList, decoyList


def read_csv(csvfile, options):
    """
	Read csv and return molList, a list of mol objects
	"""

    # open file or exit
    name, ext = os.path.splitext(csvfile)
    try:
        if ext == '.gz':
            f = gzip.open(csvfile, 'rb')
        else:
            f = open(csvfile, 'rU')
    except IOError:
        print(" \n '{f}' could not be opened\n".format(f=os.path.basename(csvfile)))
        sys.exit(1)

    # read file
    csv_reader = csv.reader(f)
    molList = []
    linenumber = 1

    for line in csv_reader:
        # get column labels from the first line
        if linenumber == 1:
            prop_indices = read_header(line, options)

        # otherwise read line & append to MolList
        else:
            mol = Molecule()
            mol = read_line(line, options, prop_indices, mol)
            # if the line's junk, skip it
            if mol == 1:
                print(" skipping molecule 'm'\n".format(m=(linenumber - 1)))
            else:
                molList.append(mol)

        linenumber += 1

    # return molList

    return molList


def read_line(csv_contents, options, prop_indices, mol):
    """
	read csv line
	"""

    # set input variables
    status_field = options.status_field
    active_label = options.active_label
    decoy_label = options.decoy_label

    # confirm the active/decoy labels have appropriate values
    # otherwise skip molecule
    active_value_matcher = re.compile(active_label)
    decoy_value_matcher = re.compile(decoy_label)
    status_label_index = prop_indices[status_field]

    if not active_value_matcher.match(csv_contents[status_label_index]) and not decoy_value_matcher.match(
            csv_contents[status_label_index]):
        print("\n molecule lacks appropriate status label")
        return 1

    # loop over property values
    for prop_label in prop_indices.keys():
        # get property value
        value_index = prop_indices[prop_label]
        prop_value = csv_contents[value_index]

        # set corresponding molecule attribute
        mol.SetProp(prop_label, prop_value)

    # return mol
    return mol


def read_header(headerlabels, options):
    """
	read csv header, returns prop_indices dictionary, {'label' : index}
	where label is the csv column label, and index the column index
	"""

    # set input variables
    status_field = options.status_field

    # confirm that status_field exists
    status_field_matcher = re.compile(status_field)
    if not any(x for x in headerlabels if status_field_matcher.match(x)):
        print("\n The status_field column '{s}' doesn't exist\n".format(s=status_field))
        sys.exit(1)

    # initiate prop_column dictionary and csv column index
    prop_indices = {}
    index = 0

    for label in headerlabels:
        prop_indices[label] = index
        index += 1

    return prop_indices


class Molecule:
    """
	csv molecule class
	"""

    def __init__(self):
        """
		class attributes
		"""
        self.id = None
        self.status = None
        self.bm = None
        self.graph = None
        self.properties = {}
        self.tags = ('id', 'status', 'bm', 'graph')

    def SetProp(self, prop, value):
        """
		set attribute
		"""
        if prop == 'id':
            self.id = value
        elif prop == 'status':
            self.status = value
        elif prop == 'bm':
            self.bm = value
        elif prop == 'graph':
            self.graph = value
        else:
            self.properties[prop] = value

    def GetProp(self, prop):
        """
		get attribute
		"""

        if prop == 'id':
            return self.id
        elif prop == 'status':
            return self.status
        elif prop == 'bm':
            return self.bm
        elif prop == 'graph':
            return self.graph
        else:
            return self.properties[prop]


class SplitInput:
    """
	Return command line input
	"""

    def __init__(self, itf):
        self.input = itf.get_string('--input')
        self.status_field = 'status'
        self.active_label = '1'
        self.decoy_label = '0'
        self.training_fraction = itf.get_float('--training_fraction')
        self.decoy_to_active = itf.get_int('--decoy_to_active')

        # set defaults for optional arguments
        if self.training_fraction is None:
            self.training_fraction = 0.5