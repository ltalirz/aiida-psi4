# -*- coding: utf-8 -*-
"""
Calculations provided by aiida_psi4.

Register calculations via the "aiida.calculations" entry point in setup.json.
"""
import json
import io
from aiida.common import datastructures, InputValidationError
from aiida.engine import CalcJob
from aiida.plugins import DataFactory
from aiida import orm

AtomicInput = DataFactory('psi4.atomic_input')


class Psi4Calculation(CalcJob):
    """
    AiiDA calculation plugin wrapping the psi4 executable.

    """
    _DEFAULT_INPUT_FILENAME = 'in.json'
    _DEFAULT_OUTPUT_FILENAME = 'out.json'

    @classmethod
    def define(cls, spec):
        """Define inputs and outputs of the calculation."""
        # yapf: disable
        super(Psi4Calculation, cls).define(spec)

        # set default values for AiiDA options
        spec.inputs['metadata']['options']['resources'].default = {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 1,
        }
        spec.inputs['metadata']['options']['parser_name'].default = 'psi4'

        # new ports
        spec.input('metadata.options.output_filename', valid_type=str, default=cls._DEFAULT_OUTPUT_FILENAME)
        spec.input('qcschema', valid_type=(orm.Dict, AtomicInput), help='Psi4 input in QCSchema JSON format')
        spec.output('qcschema', valid_type=orm.Dict, help='Psi4 output in QCSchema JSON format')

        spec.exit_code(100, 'ERROR_MISSING_OUTPUT_FILES', message='Calculation did not produce all expected output files.')
        spec.exit_code(101, 'ERROR_CALCULATION_FAILED', message='Psi4 reported calculation as unsuccessful.')



    def prepare_for_submission(self, folder):
        """
        Create input files.

        :param folder: an `aiida.common.folders.Folder` where the plugin should temporarily place all files
            needed by the calculation.
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        codeinfo = datastructures.CodeInfo()
        codeinfo.cmdline_params = [ '--qcschema', self._DEFAULT_INPUT_FILENAME, '--output', self._DEFAULT_OUTPUT_FILENAME ]
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.withmpi = self.inputs.metadata.options.withmpi

        # Add QCSchema input file
        with io.open(folder.get_abs_path(self._DEFAULT_INPUT_FILENAME), mode='w', encoding='utf-8') as fobj:
            try:
                fobj.write(json.dumps(self.inputs.qcschema.get_dict(), indent=2))
            except ValueError as exc:
                raise InputValidationError('Unable to write input file to disk') from exc

        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = [ ]
        calcinfo.retrieve_list = [self.metadata.options.output_filename]

        return calcinfo
