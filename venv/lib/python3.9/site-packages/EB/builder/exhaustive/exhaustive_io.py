__author__ = 'robswift'

import os, sys
import textwrap


def print_short_help():
    """
	Prints a short help message.
	"""
     # initiate TextWrapper class, which will handle all of the string formatting
    w = textwrap.TextWrapper()
    w.expand_tabs = False
    w.width = 115
    w.initial_indent = ' '
    w.subsequent_indent = '\t\t '
    cmd = "usage: optimize --input [(required) csv file] --outname [(required) string] " \
          "--ensemble_size [(required) integer] --fpf [(optional) float] --ncpu [(optional) integer]"
    print('')
    print(w.fill(cmd))
    print('')

def print_extended_help():
    """
	Prints an extended help message.
	"""
    #screen_hght, screen_wdth = os.popen('stty size', 'r').read().split()

    # initiate TextWrapper class, which will handle all of the string formatting
    w = textwrap.TextWrapper()
    w.expand_tabs = False
    w.width=85
    w.initial_indent = '\t'
    w.subsequent_indent = '\t  '

    print('')
    print(textwrap.fill("<exhaustive> Complete parameter list", initial_indent='  '))
    print('')

    cmd = "--input : (required) csv file of screened molecules."
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

    ensemble_size = "--ensemble_size : (required) Specifies the maximum ensemble size \
	to consider."
    print(w.fill(ensemble_size))
    print('')

    fpf = "--fpf : (optional) Rank ensembles by the number of unique active\
	molecules found above the specified value. Ensembles will be ranked by AUC if a value \
	of '1' is given. The default is to use the smallest defined false positive fraction."
    print(w.fill(fpf))
    print('')

    n_cpu = "--ncpu : (optional) Number of parallel processes to use. The default is \
	to use all available cores."
    print(w.fill(n_cpu))
    print('')


def get_interface(args):
    if len(args) < 6:
        print_short_help()
        return None
    else:
        # check the input flags
        flags = ('--input', '--outname', '--ensemble_size',
                 '--fpf', '--ncpu')
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
        self.required = ('--input', '--outname', '--ensemble_size')
        self.flags = ('--input', '--outname', '--ensemble_size', '--fpf', '--ncpu')

    def get_string(self, input_string):
        """
		Return string type user input
		"""

        if input_string in ('--input', '--outname'):

            # was the flag set?
            try:
                index = self.args.index(input_string) + 1
            except ValueError:
                # it wasn't, so if it's required, exit
                if input_string in self.required:
                    print("\n {flag} is required".format(flag=input_string))
                    print_short_help()
                    sys.exit(1)
                # it wasn't, if its optional, return the default
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
            elif input_string == '--outname':
                return format(self.args[index])

    def get_int(self, input_string):
        """
		Return integer type user input
		"""

        if input_string in ('--ensemble_size', '--ncpu'):

            # was the flag set?
            try:
                index = self.args.index(input_string) + 1
            except ValueError:
                # it wasn't, so if it's required, exit
                if input_string in self.required:
                    print("\n {flag} is required\n".format(flag=input_string))
                    print_short_help()
                    sys.exit(1)
                # otherwise return the appropriate default
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

            # a value was set, so check if its the correct type
            try:
                value = int(self.args[index])
            except ValueError:
                print("\n {flag} must be an integer".format(flag=input_string))
                print_short_help()
                sys.exit(1)

            # verify the value provided is not negative
            if value < 0:
                print("\n {flag} must be an integer greater than 0".format(flag=input_string))
                print_short_help()
                sys.exit(1)

            # everything checks out, so return the value
            return value

    def get_float(self, input_string):
        """
		Return float type user input
		"""

        if input_string == '--fpf':

            # was the flag set?
            try:
                index = self.args.index(input_string) + 1
            except ValueError:
                # it wasn't, fpf is optional, so return the appropriate default
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

            # a value was set, so check if its the correct type
            try:
                value = float(self.args[index])
            except ValueError:
                print("\n {flag} must be a float less than or equal to 1, e.g. 0.01".format(flag=input_string))
                print_short_help()
                sys.exit(1)
            if value > 1.0 or value < 0:
                print("\n {flag} must be a float less than or equal to 1, e.g. 0.01".format(flag=input_string))
                print_short_help()
                sys.exit(1)

        # everything checks out, so return the appropriate value
        return value