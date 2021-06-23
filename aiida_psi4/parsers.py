# -*- coding: utf-8 -*-
"""
Parsers provided by aiida_psi4.

Register parsers via the "aiida.parsers" entry point in setup.json.
"""
import json
import io
from aiida.engine import ExitCode
from aiida.parsers.parser import Parser
from aiida.plugins import CalculationFactory
from aiida.common import exceptions
from aiida import orm
from aiida_psi4.calculations import PSI4_FILENAMES, SinglefileData

Psi4Calculation = CalculationFactory('psi4')


class QCSchemaParser(Parser):
    """
    Parse output of Psi4 calculation in QCSchema format.
    """
    def __init__(self, node):
        """
        Initialize Parser instance

        Checks that the ProcessNode being passed was produced by a Psi4Calculation.

        :param node: ProcessNode of calculation
        :param type node: :class:`aiida.orm.ProcessNode`
        """
        super().__init__(node)
        if not issubclass(node.process_class, Psi4Calculation):
            raise exceptions.ParsingError('Can only parse Psi4Calculation')

    def parse(self, **kwargs):
        """
        Parse outputs, store results in database.

        :returns: an exit code, if parsing fails (or nothing if parsing succeeds)
        """
        if 'qcschema' in self.node.inputs:
            input_method = 'qcschema'
        if 'psiapi' in self.node.inputs:
            input_method = 'psiapi'
        output_filename = PSI4_FILENAMES[input_method]['output']

        # Check that folder content is as expected
        files_retrieved = self.retrieved.list_object_names()
        files_expected = [output_filename]
        # Note: set(A) <= set(B) checks whether A is a subset of B
        if not set(files_expected) <= set(files_retrieved):
            self.logger.error("Found files '{}', expected to find '{}'".format(
                files_retrieved, files_expected))
            return self.exit_codes.ERROR_MISSING_OUTPUT_FILES

        # add outputs
        self.logger.info("Parsing '{}'".format(output_filename))
        with self.retrieved.open(output_filename, 'rb') as handle:

            if input_method == 'psiapi':
                log_node = SinglefileData(file=handle,
                                          filename=output_filename)

            elif input_method == 'qcschema':
                output_dict = json.loads(handle.read())
                if not output_dict['success']:
                    return self.exit_codes.ERROR_CALCULATION_FAILED

                # remove stdout (don't want to store unparsed files in the database)
                log_node = SinglefileData(
                    # note: in python3.9 with AiiDA 2.0 this can be simplified to
                    # file=io.StrinIO(''.join(output_dict['stdout'])),
                    file=io.BytesIO(
                        bytes(''.join(output_dict['stdout']),
                              encoding='utf8')),
                    filename=PSI4_FILENAMES['qcschema']['output'])
                output_dict.pop('stdout')

                self.out('qcschema', orm.Dict(dict=output_dict))

            self.out('stdout', log_node)

        return ExitCode(0)
