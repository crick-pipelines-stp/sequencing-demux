"""
Tests for the samplesheet_check module
"""

import os
import sys
from pathlib import Path
import pytest
import unittest


# Alter python path
script_path = Path(__file__).parent / '../../templates'
sys.path.insert(0, str(script_path))

# Import test script
from samplesheet_check import main, check_samplesheet  # pylint: disable=E0401,C0413 # noqa: E402 # type: ignore

# Get script location
script_dir = os.path.dirname(os.path.abspath(__file__))


class TestSampleSheetProcessing(unittest.TestCase):
    """
    Test class
    """

    def test_run_valid_sheet_single_sample(self):
        """
        Test case
        """

        # Setup
        samplesheet_path = os.path.join(script_dir, '..', 'data', 'single_sample.csv')
        output_path = 'samplesheet.checked.csv'
        args = [
            '--process_name', 'SAMPLESHEET_CHECK',
            '--sample', samplesheet_path,
            '--output', output_path
        ]

        # Run
        main(args)

        # Assert
        self.assertTrue(os.path.exists(output_path), "Output file was not created")

        # Teardown
        if os.path.exists(output_path):
            os.remove(output_path)

    def test_error_too_many_columns(self):
        """
        Test case
        """

        samplesheet_path = os.path.join(script_dir, '..', 'data', 'error', 'too_many_columns.csv')

        with pytest.raises(SystemExit):
            check_samplesheet(samplesheet_path, "samplesheet.checked.csv")
