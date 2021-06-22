#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run a test calculation on localhost.

Usage: ./example_01.py
"""
from os import path
import click
from aiida import cmdline, engine
from aiida.plugins import DataFactory, CalculationFactory
from aiida_psi4 import helpers

INPUT_DIR = path.join(path.dirname(path.realpath(__file__)), 'input_files')


def test_run(psi4_code):
    """Run a calculation on the localhost computer.

    Uses test helpers to create AiiDA Code on the fly.
    """
    if not psi4_code:
        # get code
        computer = helpers.get_computer()
        psi4_code = helpers.get_code(entry_point='psi4', computer=computer)

    # Prepare input parameters
    AtomicInput = DataFactory('psi4')
    parameters = AtomicInput({'ignore-case': True})

    SinglefileData = DataFactory('singlefile')
    file1 = SinglefileData(file=path.join(INPUT_DIR, 'file1.txt'))
    file2 = SinglefileData(file=path.join(INPUT_DIR, 'file2.txt'))

    # set up calculation
    inputs = {
        'code': psi4_code,
        'parameters': parameters,
        'file1': file1,
        'file2': file2,
        'metadata': {
            'description': 'Test job submission with the aiida_psi4 plugin',
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # from aiida.engine import submit
    # future = submit(CalculationFactory('psi4'), **inputs)
    result = engine.run(CalculationFactory('psi4'), **inputs)

    computed_diff = result['psi4'].get_content()
    print('Computed diff between files: \n{}'.format(computed_diff))


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
def cli(code):
    """Run example.

    Example usage: $ ./example_01.py --code diff@localhost

    Alternative (creates diff@localhost-test code): $ ./example_01.py

    Help: $ ./example_01.py --help
    """
    test_run(code)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter