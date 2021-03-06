[![Build Status](https://github.com/aiidateam/aiida-psi4/workflows/ci/badge.svg?branch=master)](https://github.com/aiidateam/aiida-psi4/actions)
[![Coverage Status](https://coveralls.io/repos/github/aiidateam/aiida-psi4/badge.svg?branch=master)](https://coveralls.io/github/aiidateam/aiida-psi4?branch=master)
[![Docs status](https://readthedocs.org/projects/aiida-psi4/badge)](http://aiida-psi4.readthedocs.io/)
[![PyPI version](https://badge.fury.io/py/aiida-psi4.svg)](https://badge.fury.io/py/aiida-psi4)

# aiida-psi4

AiiDA plugin for the Psi4 Quantum Chemistry package.

Note: This plugin uses the QCSchema input format of Psi4 and therefore requires Psi4 version 1.4 or above.

## Features

 * Describe input using `AtomicInput` dictionaries and obtain pre-parsed QCSchema output
   (see [`examples/example_01.py`](examples/example_01.py))
   ```python
   AtomicInput = DataFactory('psi4.atomic_input')
   inputs['qcschema'] = AtomicInput({...})
   ...
   qcdict = result['qcschema'].get_dict()
   ```
* Provide raw python as string and obtain log file as output
  (see [`examples/example_02.py`](examples/example_02.py))
   ```python
   psiapi_str = """import psi4; psi4.geometry('''...'''); psi4.energy*"ccsd(t)/cc-pvtz")"""
   inputs['psiapi'] = Str(psiapi_str)
   ...
   log = result['stdout'].get_content()
   ```

## Installation

```shell
pip install aiida-psi4
verdi quicksetup  # better to set up a new profile
verdi plugin list aiida.calculations  # should now show your calclulation plugins
```


## Usage

Here goes a complete example of how to submit a test calculation using this plugin.

A quick demo of how to submit a calculation:
```shell
verdi daemon start     # make sure the daemon is running
cd examples
./example_01.py        # run test calculation
verdi process list -a  # check record of calculation
```

## Development

```shell
git clone https://github.com/aiidateam/aiida-psi4 .
cd aiida-psi4
pip install -e .[pre-commit,testing]  # install extra dependencies
pre-commit install  # install pre-commit hooks
pytest -v  # discover and run all tests
```

See the [developer guide](http://aiida-psi4.readthedocs.io/en/latest/developer_guide/index.html) for more information.

## License

MIT
## Contact

leopold.talirz@gmail.com
