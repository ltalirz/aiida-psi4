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
SinglefileData = DataFactory('singlefile')

PSI4_FILENAMES = {
    'qcschema': {
        'input': 'input.json',
        'output': 'output.json',
    },
    'psiapi': {
        'input': 'input.py',
        'output': 'output.log',
    }
}


class Psi4Calculation(CalcJob):
    """
    AiiDA calculation plugin wrapping the psi4 executable.

    """
    @classmethod
    def define(cls, spec):
        """Define inputs and outputs of the calculation."""
        # yapf: disable
        super().define(spec)

        # set default values for AiiDA options
        spec.inputs['metadata']['options']['resources'].default = {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 1,
        }
        spec.inputs['metadata']['options']['parser_name'].default = 'psi4'

        # new ports
        #spec.input('metadata.options.output_filename', valid_type=str, default=cls._DEFAULT_OUTPUT_FILENAME)
        spec.input('qcschema', valid_type=(orm.Dict, AtomicInput), required=False, help='Psi4 input in QCSchema JSON format')
        spec.input('psiapi', valid_type=(orm.Str,  SinglefileData),
            required=False, help='Psi4 input in PsiAPI python format')
        spec.output('qcschema', valid_type=orm.Dict, help='Psi4 output in QCSchema JSON format')
        spec.output('stdout', valid_type=SinglefileData, help='Psi4 logfile')

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
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.withmpi = self.inputs.metadata.options.withmpi


        # Add input file
        if 'qcschema' in self.inputs and 'psiapi' in self.inputs:
            raise InputValidationError('Do not mix "qcschema" and "psithon" input')

        if 'qcschema' in self.inputs:
            filenames = PSI4_FILENAMES['qcschema']
            codeinfo.cmdline_params = [ '--qcschema', filenames['input'], '--output', filenames['output']]
            input_string = json.dumps(self.inputs.qcschema.get_dict(), indent=2)
        elif 'psiapi' in self.inputs:
            filenames = PSI4_FILENAMES['psiapi']
            codeinfo.cmdline_params = [ filenames['input'], '--output', filenames['output']]
            if isinstance(self.inputs.psiapi, SinglefileData):
                with self.inputs.psiapi.open('r') as handle:
                    input_string = handle.read()
            elif isinstance(self.inputs.psiapi, orm.Str):
                input_string = self.inputs.psiapi.value
        else:
            raise InputValidationError('Please provide either "qcschema" or "psithon" input')

        # write input file
        with io.open(folder.get_abs_path(filenames['input']), mode='w', encoding='utf-8') as fobj:
            try:
                fobj.write(input_string)
            except ValueError as exc:
                raise InputValidationError(f'Unable to write input file {filenames["input"]} to disk') from exc

        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = [ ]
        calcinfo.retrieve_list = [filenames['output']]

        return calcinfo
