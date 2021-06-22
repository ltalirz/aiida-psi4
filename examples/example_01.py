#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run a test calculation on localhost.

Usage: ./example_01.py
"""
from os import path
from pprint import pprint
import click
import pytest
from aiida import cmdline, engine
from aiida.plugins import DataFactory, CalculationFactory

INPUT_DIR = path.join(path.dirname(path.realpath(__file__)), 'input_files')


def test_run(psi4_code):
    """Run a calculation on the localhost computer.

    Uses test helpers to create AiiDA Code on the fly.
    """
    # Prepare input parameters
    TEST_DICT = {
        'schema_name': 'qcschema_input',
        'schema_version': 1,
        'molecule': {
            'geometry': [
                0.0, 0.0, -0.1294769411935893, 0.0, -1.494187339479985,
                1.0274465079245698, 0.0, 1.494187339479985, 1.0274465079245698
            ],
            'symbols': ['O', 'H', 'H']
        },
        'driver': 'energy',
        'model': {
            'method': 'CCSD(T)',
            'basis': '6-31g'
        },
        'keywords': {
            'scf_type': 'df',
            'mp2_type': 'df',
            'cc_type': 'df',
            'scf_properties': ['mayer_indices']
        }
    }
    AtomicInput = DataFactory('psi4.atomic_input')
    atomic_input = AtomicInput(TEST_DICT)

    # set up calculation
    inputs = {
        'code': psi4_code,
        'qcschema': atomic_input,
        'metadata': {
            'description': 'Test job submission with the aiida_psi4 plugin',
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # from aiida.engine import submit
    # future = submit(CalculationFactory('psi4'), **inputs)
    result = engine.run(CalculationFactory('psi4'), **inputs)
    qcdict = result['qcschema'].get_dict()
    pprint(qcdict)
    assert qcdict['success']
    assert qcdict['return_energy'] == pytest.approx(-76.120695408714, 0.0001)


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
