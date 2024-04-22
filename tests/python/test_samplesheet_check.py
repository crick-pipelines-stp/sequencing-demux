"""
Tests for the samplesheet_check module
"""

import sys
from pathlib import Path
import pytest


# Alter python path
script_path = Path(__file__).parent / '../../modules/local/samplesheet_check/templates'
sys.path.insert(0, str(script_path))

# Import test script
from samplesheet_check import main, check_samplesheet  # pylint: disable=E0401,C0413 # noqa: E402 # type: ignore


def test_run_valid_sheet_single_sample():
    """
    Test case
    """

    samplesheet_path = "tests/data/samplesheets/single_sample.csv"
    args = [
        '--process_name', 'SAMPLESHEET_CHECK',
        '--sample', samplesheet_path,
        '--output', 'samplesheet.checked.csv'
    ]

    main(args)


def test_error_too_many_columns():
    """
    Test case
    """

    samplesheet_path = "tests/data/samplesheets/error/too_many_columns.csv"

    with pytest.raises(SystemExit):
        check_samplesheet(samplesheet_path, "samplesheet.checked.csv")
