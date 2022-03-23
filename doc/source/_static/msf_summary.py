#! /bin/env python

import sys
import os
import argparse
import pandas as pd



def summary_best_solution(best_solution):
    """

    :param best_solution: the data contained in best_solution file
    :type best_solution: :class:`pandas.DataFrame` object
    :return: the summary of the best_solution data
    :rtype: :class:`pandas.DataFrame` object
    """
    selection = best_solution[['replicon', 'sys_id', 'model_fqn']]
    dropped = selection.drop_duplicates(subset=['replicon', 'sys_id'])
    summary = pd.crosstab(index=dropped.replicon, columns=dropped['model_fqn'])
    return summary


def summary_all_best_solution(all_best_solution):
    """

    :param all_best_solution: the data contained in all_best_solutions file
    :type best_solution: :class:`pandas.DataFrame` object
    :return: the summary of the all_best_solutions data
    :rtype: :class:`pandas.DataFrame` object
    """
    selection = all_best_solution[['sol_id', 'replicon', 'sys_id', 'model_fqn']]
    dropped = selection.drop_duplicates(subset=['sol_id', 'replicon', 'sys_id'])
    summary = pd.crosstab(index=[dropped.sol_id, dropped.replicon], columns=dropped['model_fqn'])
    return summary


def parse_args(args):
    """
    Pars arguments on command line

    :param args: the command lines arguments
    :type args: list of strings
    :return: the parsed comandline arguments
    :rtype: :class:`argparse.Namespace` object
    """
    parser = argparse.ArgumentParser(
        description="compute a quick summary for best_solution.tsv or all_best_solutions.tsv file",
        epilog="For more details, visit the MacSyFinder website and see the MacSyFinder documentation.")
    parser.add_argument('best_solution',
                        help="The path to a best_solution.tsv or all_best_solutions.tsv")
    parsed_args = parser.parse_args(args)
    return parsed_args


def main():
    """
    main entry point, read the data, send them to the right function
    depending if it's a best_solution ot all_best_solution data
    writing the results on file.
    """
    args = parse_args(sys.argv[1:])
    if not os.path.exists(args.best_solution):
        raise FileNotFoundError(f"'{args.best_solution}': No such file.")
    solution = pd.read_csv(args.best_solution, sep='\t', comment='#')

    if 'sol_id' in solution.columns:
        summary = summary_all_best_solution(solution)
    else:
        summary = summary_best_solution(solution)

    base , ext = os.path.splitext(args.best_solution)
    out = f"{base}-summary{ext}"
    summary.to_csv(out, sep='\t')
    print(f"summary is save in '{out}'")


if __name__ == '__main__':
    main()