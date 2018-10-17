import unittest
import sys
import os


def _run(test_files, test_root_path, verbosity=0):
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


def run_unittests(test_files, verbosity=0):
    """
    Execute Unit Tests

    :param test_files: the file names of tests to run.
    of it is empty, discover recursively tests from 'tests/unit' directory.
    a test is python module with the test_*.py pattern
    :type test_files: list of string
    :param verbosity: the verbosity of the output
    :type verbosity: positive int
    :return: True if the test passed successfully, False otherwise.
    :rtype: bool
    """
    test_root_path = os.path.join('tests', 'unit')
    return _run(test_files, test_root_path, verbosity)


def run_integration_tests(test_files, verbosity=0):
    """
    Execute Integration Tests

    :param test_files: the file names of tests to run.
    of it is empty, discover recursively tests from 'tests/unit' directory.
    a test is python module with the test_*.py pattern
    :type test_files: list of string
    :param verbosity: the verbosity of the output
    :type verbosity: positive int
    :return: True if the test passed successfully, False otherwise.
    :rtype: bool
    """
    test_root_path = os.path.join('tests', 'integration')
    return _run(test_files, test_root_path, verbosity)


def run_functional_tests(test_files, verbosity=0):
    """
    Execute Functional Tests

    :param test_files: the file names of tests to run.
    of it is empty, discover recursively tests from 'tests/unit' directory.
    a test is python module with the test_*.py pattern
    :type test_files: list of string
    :param verbosity: the verbosity of the output
    :type verbosity: positive int
    :return: True if the test passed successfully, False otherwise.
    :rtype: bool
    """
    test_root_path = os.path.join('tests', 'functional')
    return _run(test_files, test_root_path, verbosity)


if __name__ == '__main__':



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
        print("\n", "#" * 70, sep='')
        print("Test Runner: Unit tests")
        print("#" * 70)

        old_path = sys.path
        if 'MACSY_HOME' in os.environ:
            MACSY_HOME = os.environ['MACSY_HOME']
            if MACSY_HOME not in sys.path:
                sys.path.insert(0, MACSY_HOME)
        unit_results = run_unittests(args.tests, verbosity=args.verbosity)
        result_all_tests.append(unit_results)
        sys.path = old_path

    if args.functional:
        print("\n", "#" * 70, sep='')
        print("Test Runner: Functional tests")
        print("#" * 70)

        old_path = sys.path
        if 'MACSY_HOME' in os.environ:
            MACSY_HOME = os.environ['MACSY_HOME']
            if MACSY_HOME not in sys.path:
                sys.path.insert(0, MACSY_HOME)
        else:
            home_tests = os.path.abspath(os.path.join(os.path.dirname(__file__)))
            # we are in the case where we tests an installed macsyfinder
            # so the libraries are already in PYTHONPATH
            # but test are not
            # we must had tests parent dir in pythonpath
            # but after the standard libraries containing macsyfinder
            # as we want to run macsyfinder using installed libraries
            sys.path.append(home_tests)
        functional_results = run_functional_tests(args.tests, verbosity=args.verbosity)
        result_all_tests.append(functional_results)
        sys.path = old_path

    if args.integration:
        print("\n", "#" * 70, sep='')
        print("Test Runner: Integration tests")
        print("#" * 70)
        old_path = sys.path
        integration_results = run_integration_tests(args.tests, verbosity=args.verbosity)
        result_all_tests.append(integration_results)
        sys.path = old_path

    if all(result_all_tests):
        sys.exit(0)
    else:
        sys.exit(1)
