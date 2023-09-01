__author__ = 'robswift'


import os, sys
import textwrap


def print_short_help():
    """
    print a short help message
    :return:
    """
    w = textwrap.TextWrapper()
    w.expand_tabs = False
    w.width = 125
    w.initial_indent = ' '
    w.subsequent_indent = '\t\t     '

    input = "usage: postanalysis --input [(required) csv] --outname [(required) string]\
	--ensemble_list [(required) 1st.csv ... Nth.csv] --compare [(optional)]\
	--fpf [(optional) float] --plot [(optional)] --write_roc [(optional)]"
    print('')
    print(w.fill(input))
    print('')


def print_extended_help():
    """
    print a detailed help message
    :return:
    """
    w = textwrap.TextWrapper()
    w.expand_tabs = False
    w.width = 85
    w.initial_indent = '\t'
    w.subsequent_indent = '\t  '

    print('')
    print(textwrap.fill("<postanalyze> Complete parameter list:", initial_indent=''))
    print('')

    cmd = "--input : (required) csv file to split into training and test sets"
    print(w.fill(cmd))
    cmd = "\t\tColumns should be as follows:"
    print(w.fill(cmd))
    print('')
    cmd="\t\t    id, status, receptor_1, receptor_2, ..., receptor_N"
    print(w.fill(cmd))
    cmd="\t\t  CH44,      1,       -9.7,       -9.3, ...,      -10.2"
    print(w.fill(cmd))
    cmd="\t\t  ZN44,      0,       -6.6,       -6.1, ...,       -6.8"
    print(w.fill(cmd))
    print('')
    cmd="\t\tid is a unique molecular identifier"
    print(w.fill(cmd))
    cmd="\t\tstatus takes a value of '1' if the molecule is active and '0' otherwise."
    print(w.fill(cmd))
    cmd="\t\treceptor_1 through receptor_N are docking scores."
    print(w.fill(cmd))
    print('')

    outname = "--outname : (required) the prefix of the outputfiles."
    print(w.fill(outname))
    print('')

    ensemble_list = "--ensemble_list : (required) a list of csv files that contain the queries in\
	the ensemble. For example, 'Ensemble_1_queries.csv Ensemble_2_queries.csv\
	Ensemble_3_queries.csv ...'"
    print(w.fill(ensemble_list))
    print('')

    compare = "--compare : (optional) Compare the virtual screening results for the\
	ensembles specified after the '--ensemble_list' flag. No more than two ensembles\
	may be specified at once."
    print(w.fill(compare))
    print('')

    fpf = "--fpf : (optional) Evaluate ensemble performance at the set of specified FPF values. \
    By default, values of '0.0001', '0.001', '0.01', and '0.05' are considered, if they are defined."
    print(w.fill(fpf))
    print('')

    plot = "--plot : (optional) Generate ROC plots of the input ensembles and their\
	members."
    print(w.fill(plot))
    print('')

    roc_data = "--write_roc : (optional) if the '--write_roc' flag is set, a 'ROC_DATA' \
	directory will be created, & ROC data points will be written there for each\
	ensemble. The default is not to write ROC data points."
    print(w.fill(roc_data))
    print('')


def get_interface(args):
    if len(args) < 4:
        print_short_help()
        return None
    else:
        # check the input flags
        flags = ('--input', '--ensemble_list', '--fpf', '--outname',
                 '--plot', '--write_roc', '--compare')
        for inputarg in [x for x in args if '--' in x]:
            if inputarg not in flags:
                print("\n Unrecognized input flag: {flag}\n".format(flag=inputarg))
                print_short_help()
                sys.exit(1)
        parser = ParseArgs(args)
        return parser


class ParseArgs:
    def __init__(self, args):
        self.args = args
        self.required = ('--input', '--ensemble_list', '--outname')
        self.flags = ('--input', '--ensemble_list', '--write_roc', '--outname', '--fpf',
                      '--plot', '--compare')

    def get_list(self, input_string):
        """
        Return a list of user input
        :param input_string:
        :return:
        """

        if input_string in ('--ensemble_list', '--fpf'):

            # was the flag set?
            try:
                index_low = self.args.index(input_string) + 1
            except ValueError:
                if input_string in self.required:
                    print("\n {flag} is required".format(flag=input_string))
                    print_short_help()
                    sys.exit(1)
                else:
                    return None

            # the flag was set, so check if a value was set, otherwise exit
            try:
                if self.args[index_low] in self.flags:
                    print("\n {flag} was set but a value was not specified".format(flag=input_string))
                    print_short_help()
                    sys.exit(1)
            except IndexError:
                print("\n {flag} was set but a value was not specified".format(input_string))
                print_short_help()
                sys.exit(1)

            # at least one value was set
            index_high = index_low
            try:
                # if the flag wasn't the last argument specified
                while self.args[index_high] not in self.flags:
                    index_high += 1
            except IndexError:
                # if it was, then handle it accordingly
                index_high = self.args.index(self.args[-1])
                return self.args[index_low:index_high + 1]
            # return a list of input files
            if index_low == index_high:
                inputList = []
                inputList.append(self.args[index_low])
                return inputList
            else:
                return self.args[index_low:index_high]

    def get_string(self, input_string):
        """
		Return string type user input
		"""

        if input_string in ('--input', '--outname'):

            # was the flag set?
            try:
                index = self.args.index(input_string) + 1
            except ValueError:
                # was it required? If so, print error and exit
                if input_string in self.required:
                    print("\n {flag} is required".format(flag=input_string))
                    print_short_help()
                    sys.exit(1)
                # if not required, set and return default value
                else:
                    return None

            # the flag was set, so check if a value was set, otherwise exit
            try:
                if self.args[index] in self.flags:
                    print("\n {flag} was set but a value was not specified".format(flag=input_string))
                    print_short_help()
                    sys.exit(1)
            except IndexError:
                print("\n {flag} was set but a value was not specified".format(flag=input_string))
                print_short_help()
                sys.exit(1)

            # a value was set, so check and assign the appropriate value or exit
            if input_string == '--input':
                return os.path.abspath(self.args[index])
            else:
                return self.args[index]

    def get_boolean(self, input_string):
        """
		Return boolean type user input
		"""

        if input_string in ('--write_roc', '--plot', '--compare'):

            # was the flag set?
            try:
                index = self.args.index(input_string) + 1
            except ValueError:
                # it wasn't, args are optional, so return the appropriate default
                return False

            # the flag was set, so return the True
            return True