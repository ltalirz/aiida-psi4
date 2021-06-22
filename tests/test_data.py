# -*- coding: utf-8 -*-
""" Tests for data

"""
import pytest
from pydantic.error_wrappers import ValidationError  # pylint: disable=no-name-in-module
from aiida.plugins import DataFactory
# from . import TEST_DIR

AtomicInput = DataFactory('psi4.atomic_input')

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


def test_atomic_input():
    """Test the AtomicInput Data class
    """
    # Passing valid dict should not raise
    AtomicInput(TEST_DICT)

    # Adding an invalid key should raise pydantic validation error
    t2 = TEST_DICT.copy()
    t2['invalid_key'] = 123
    with pytest.raises(ValidationError):
        AtomicInput(t2)
