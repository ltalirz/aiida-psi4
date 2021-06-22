# -*- coding: utf-8 -*-
"""pytest fixtures for simplified testing."""
from __future__ import absolute_import
import pytest
from tests import DATA_DIR

pytest_plugins = [
    'aiida.manage.tests.pytest_fixtures', 'aiida_testing.mock_code'
]  # pylint: disable=invalid-name


@pytest.fixture(scope='function', autouse=True)
def clear_database_auto(clear_database):  # pylint: disable=unused-argument
    """Automatically clear database in between tests."""


@pytest.fixture(scope='function')
def psi4_code(mock_code_factory):
    """Create mocked "psi4" code."""
    return mock_code_factory(
        label='psi4-1.4rc2',
        data_dir_abspath=DATA_DIR,
        entry_point='psi4',
        # files *not* to copy into the data directory
        ignore_paths=('_aiidasubmit.sh', ))
