__author__ = 'robswift'

import sys

app_names = ["exhaustive", "fastheuristic", "slowheuristic",  "postanalysis", "splitter"]

apps = {
    "exhaustive" : ["exhaustive.exhaustive_io", "exhaustive.exhaustive",
                 "\t  Determine the best performer by considering all possible ensembles, O(2^N)."],
    "fastheuristic" : ["fastheuristic.fastheuristic_io", "fastheuristic.fastheuristic",
                    "Determine the best ensemble using a O(N) heuristic."],
    "slowheuristic" : ["slowheuristic.slowheuristic_io", "slowheuristic.slowheuristic",
                    "Determine the best ensemble using an O(N^2) heuristic."],
    "postanalysis" : ["postanalysis.postanalysis_io", "postanalysis.postanalysis",
                    "Plot and analyze the performance of one or more ensembles."],
    "splitter" : ["splitter.splitter_io", "splitter.splitter",
                  "\t  Split csv input into training and test sets."],
}

help_string = """
usage:	ensemblebuilder <mode> <args>

	ensemblebuilder <mode> to show help for that mode.

modes: """


def main(argv=[__name__]):
    try:
        import numpy
    except ImportError:
        print("NumPy must be installed")
        print("try running: pip install numpy")
        return 1
    try:
        import scipy
    except ImportError:
        print("SciPy must be installed")
        print("try: pip install scipy")
        return 1
    if len(argv) == 1 or argv[1] not in apps:
        print(help_string)
        for mode in app_names:
            print("    %s\t  %s" % (mode, apps[mode][2]))
        print
        return 1

    mode = argv[1]
    io_stream_name = 'EB.builder.' + apps[mode][0]
    builder_name = 'EB.builder.' + apps[mode][1]

    __import__(io_stream_name)
    __import__(builder_name)

    io_module = sys.modules[io_stream_name]
    builder_module = sys.modules[builder_name]

    if len(argv) == 2:
        io_module.print_extended_help()
        return 1
    else:
        run_args = argv[2:]

    itf = io_module.get_interface(run_args)

    return builder_module.run(itf)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
