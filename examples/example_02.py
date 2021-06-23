#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run a test calculation using PsiAPI input.

Usage: ./example_02.py --code psi4
"""
from pprint import pprint
import click
from aiida import cmdline, engine, orm
from aiida.plugins import CalculationFactory


def test_run(psi4_code):
    """Run a calculation using PsiAPI input.
    """
    # Prepare input parameters
    psiapi_str = '''
import psi4
psi4.geometry("""
Ne
Ne 1 3.0
""")
psi4.set_options({"freeze_core": "True"})
psi4.energy("ccsd(t)/cc-pvtz")
'''

    # set up calculation
    inputs = {
        'code': psi4_code,
        'psiapi': orm.Str(psiapi_str),
        'metadata': {
            'description': 'Test job submission with the aiida_psi4 plugin',
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # from aiida.engine import submit
    # future = submit(CalculationFactory('psi4'), **inputs)
    result = engine.run(CalculationFactory('psi4'), **inputs)
    log = result['stdout'].get_content()
    pprint(log)
    # CCSD total energy
    assert '-257.59' in log


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
def cli(code):
    """Run example.

    Example usage: $ ./example_02.py --code psi4@localhost

    Help: $ ./example_02.py --help
    """
    test_run(code)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
