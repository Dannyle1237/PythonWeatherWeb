__author__ = 'robswift'

import os
import sys
import textwrap


def print_short_help():
    """
	Prints a short help message.
	"""
    # initiate TextWrapper class, which will handle all of the string formatting
    w = textwrap.TextWrapper()
    w.expand_tabs = False
    w.width = 130
    w.initial_indent = ' '
    w.subsequent_indent = '\t\t'

    input = "usage: split --input [(required) csv file] --training_fraction [(optional) fraction]\
	--decoy_to_active [(optional) integer]"
    print('')
    print(w.fill(input))
    print('')


def print_extended_help():
    """
	Prints an extended help message.
	"""
    # initiate TextWrapper class, which will handle all of the string formatting
    w = textwrap.TextWrapper()
    w.expand_tabs = False
    w.width=110
    w.initial_indent = '    '
    w.subsequent_indent = '      '

    print('')
    print(textwrap.fill("<split> Complete parameter list:", initial_indent=''))
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

    tfrac = "--training_fraction : (optional) The fraction of input active molecules\
	allocated to the training set, e.g. 0.50. Defaults to allocate half to the training\
	set."
    print(w.fill(tfrac))
    print('')

    d2a = "--decoy_to_active : (optional) The decoy to active ratio to establish in the \
	training and validation sets. Defaults to maintain the input file ratio."
    print(w.fill(d2a))
    print('')


def get_interface(args):
    if len(args) < 2:
        print_short_help()
        return None
    else:
        # check the input flags
        flags = ('--input', '--training_fraction', '--decoy_to_active',)
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
        self.required = ('--inputfile',)
        self.flags = ('--input', '--training_fraction', '--decoy_to_active',)

    def get_string(self, input_string):
        """
		Return string type user input
		"""

        if input_string in ('--input',):

            # was the flag set?
            try:
                index = self.args.index(input_string) + 1
            except ValueError:
                # it wasn't set, so if it's required, exit
                if input_string in self.required:
                    print("\n {flag} is required".format(flag=input_string))
                    print_short_help()
                    sys.exit(1)
                # otherwise return the default
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

            # a value was set, so return the appropriate value
            return format(self.args[index])

    def get_float(self, input_string):
        """
		Return float type user input
		"""

        if input_string == '--training_fraction':

            # was the flag set?
            try:
                index = self.args.index(input_string) + 1
            except ValueError:
                # it wasn't, it's optional, so return the appropriate default
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
                print("\n {flag} must be a float less than or equal to 1, e.g. 0.4".format(flag=input_string))
                print_short_help()
                sys.exit(1)
            if value > 1.0 or value < 0:
                print("\n {flag} must be a float less than or equal to 1, e.g. 0.4".format(flag=input_string))
                print_short_help()
                sys.exit(1)

            # everything checks out, so return the appropriate value
            return value

    def get_int(self, input_string):
        """
		Return integer type user input
		"""

        if input_string in ('--decoy_to_active'):

            # was the flag set?
            try:
                index = self.args.index(input_string) + 1
            except ValueError:
                # it wasn't, so return the default
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

            # verify the value provided is not less than 1
            if value < 1:
                print("\n {flag} must be an integer greater than or equal to 1".format(flag=input_string))
                print_short_help()
                sys.exit(1)

            # everything checks out, so return the value
            return value
