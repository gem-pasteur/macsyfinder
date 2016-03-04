import unittest
import sys
import os


def _run(lib, test_files, test_root_path, verbosity=0):
    if lib not in sys.path:
        sys.path.insert(0, lib)
    if not test_files:
        suite = unittest.TestLoader().discover(test_root_path, pattern="test_*.py")
    else:
        test_files = [t for t in test_files if test_root_path in t]
        suite = unittest.TestSuite()
        for test_file in test_files:
            if os.path.exists(test_file):
                if os.path.isfile(test_file):
                    fpath, fname = os.path.split(test_file)
                    suite.addTests(unittest.TestLoader().discover(fpath, pattern=fname))
                elif os.path.isdir(test_file):
                    suite.addTests(unittest.TestLoader().discover(test_file))
            else:
                sys.stderr.write("{0} : no such file or directory\n".format(test_file))

    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)
    return result.wasSuccessful()


def run_unittests(lib, test_files, verbosity=0):
    """
    Execute Unit Tests

    :param lib: the path where is the macsypy
    :type lib: string
    :param test_files: the file names of tests to run.
    of it is empty, discover recursively tests from 'tests/unit' directory.
    a test is python module with the test_*.py pattern
    :type test_files: list of string
    :param verbosity: the verbosity of the output
    :type verbosity: positive int
    :return: True if the test passed successfully, False otherwise.
    :rtype: bool
    """
    test_root_path = 'tests/unit'
    return _run(lib, test_files, test_root_path, verbosity)


def run_integration_tests(lib, test_files, verbosity=0):
    """
    Execute Integration Tests

    :param lib: the path where is the macsypy
    :type lib: string
    :param test_files: the file names of tests to run.
    of it is empty, discover recursively tests from 'tests/unit' directory.
    a test is python module with the test_*.py pattern
    :type test_files: list of string
    :param verbosity: the verbosity of the output
    :type verbosity: positive int
    :return: True if the test passed successfully, False otherwise.
    :rtype: bool
    """
    test_root_path = 'tests/integration'
    return _run(lib, test_files, test_root_path, verbosity)


def run_functional_tests(lib, test_files, verbosity=0):
    """
    Execute Functional Tests

    :param lib: the path where is the macsypy
    :type lib: string
    :param test_files: the file names of tests to run.
    of it is empty, discover recursively tests from 'tests/unit' directory.
    a test is python module with the test_*.py pattern
    :type test_files: list of string
    :param verbosity: the verbosity of the output
    :type verbosity: positive int
    :return: True if the test passed successfully, False otherwise.
    :rtype: bool
    """
    test_root_path = 'tests/functional'
    return _run(lib, test_files, test_root_path, verbosity)


if __name__ == '__main__':
    if 'MACSY_HOME' in os.environ:
        MACSY_HOME = os.environ['MACSY_HOME']
    else:
        MACSY_HOME = os.path.abspath(os.path.join(os.path.dirname(__file__)))

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("tests",
                        nargs='*',
                        default=False,
                        help="name of test to execute")

    parser.add_argument("--unit",
                        dest='unit',
                        action='store_true',
                        default=False,
                        help="execute unit tests")

    parser.add_argument("--functional",
                        dest='functional',
                        action='store_true',
                        default=False,
                        help="execute functional tests")

    parser.add_argument("--integration",
                        dest='integration',
                        action='store_true',
                        default=False,
                        help="execute integration tests")

    parser.add_argument("-v", "--verbose",
                        dest="verbosity",
                        action="count",
                        help="set the verbosity level of output",
                        default=0
                        )

    args = parser.parse_args()
    if not any((args.unit, args.integration, args.functional)):
        args.unit, args.integration, args.functional = True, True, True

    result_all_tests = []

    if args.unit:
        print "#" * 70
        print "Test Runner: Unit tests"
        print "#" * 70
        unit_results = run_unittests(MACSY_HOME, args.tests, verbosity=args.verbosity)
        result_all_tests.append(unit_results)

    if args.functional:
        print "#" * 70
        print "Test Runner: Functional tests"
        print "#" * 70
        functional_results = run_functional_tests(MACSY_HOME, args.tests, verbosity=args.verbosity)
        result_all_tests.append(functional_results)

    if args.integration:
        print "#" * 70
        print "Test Runner: Integration tests"
        print "#" * 70
        integration_results = run_integration_tests(MACSY_HOME, args.tests, verbosity=args.verbosity)
        result_all_tests.append(integration_results)

    if all(result_all_tests):
        sys.exit(0)
    else:
        sys.exit(1)