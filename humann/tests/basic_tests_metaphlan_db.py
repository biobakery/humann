"""
Tests for MetaPhlAn database selection and verification in humann.search.prescreen
"""

import unittest
import logging
import os
import shutil
import tempfile

import cfg
import utils

from humann import config
from humann.search import prescreen


class TestMetaphlanDatabaseCheck(unittest.TestCase):
    """
    Test prescreen.verify_metaphlan_db_version
    """

    def setUp(self):
        logging.getLogger("humann.search.prescreen").addHandler(logging.NullHandler())
        self.tempdir = tempfile.mkdtemp(prefix="humann_mpa_stub_")

    def tearDown(self):
        shutil.rmtree(self.tempdir, ignore_errors=True)

    def stub_metaphlan(self, version_output):
        """Write an executable standing in for 'metaphlan --version'"""

        path = os.path.join(self.tempdir, "metaphlan_stub")
        with open(path, "w") as handle:
            handle.write("#!/bin/sh\ncat <<'EOF'\n" + version_output + "\nEOF\n")
        os.chmod(path, 0o755)
        return path

    def test_accepts_the_required_database(self):
        """The required database is found among the installed databases"""

        exe = self.stub_metaphlan(
            "MetaPhlAn version 4.2.6 (28 Jul 2026)\n"
            "Installed databases: mpa_vOct22_CHOCOPhlAnSGB_202403")
        prescreen.verify_metaphlan_db_version("vOct22_CHOCOPhlAnSGB_202403", exe)

    def test_rejects_a_different_database(self):
        """A MetaPhlAn reporting only other databases is rejected"""

        exe = self.stub_metaphlan(
            "MetaPhlAn version 4.2.6 (28 Jul 2026)\n"
            "Installed databases: mpa_vJan25_CHOCOPhlAnSGB_202503")
        with self.assertRaises(SystemExit):
            prescreen.verify_metaphlan_db_version("vOct22_CHOCOPhlAnSGB_202403", exe)

    def test_accepts_metaphlan_that_does_not_report_databases(self):
        """MetaPhlAn releases that print no database line are not rejected

        Only 4.1.2 and 4.2.4+ list installed databases. On 4.0.x and 4.1.1
        --version prints the version alone, and treating that as a mismatch
        rejected every supported MetaPhlAn except 4.1.2.
        """

        exe = self.stub_metaphlan("MetaPhlAn version 4.1.1 (23 May 2024)")
        prescreen.verify_metaphlan_db_version("vOct22_CHOCOPhlAnSGB_202403", exe)

    def test_database_line_is_recognised(self):
        """The installed-databases line is picked out of the version output"""

        lines = prescreen.db_lines(
            "MetaPhlAn version 4.2.6 (28 Jul 2026)\n"
            "Installed databases: mpa_vOct22_CHOCOPhlAnSGB_202403")
        self.assertEqual(len(lines), 1)
        self.assertTrue(lines[0].startswith("Installed databases:"))

    def test_missing_metaphlan_is_an_error(self):
        """A MetaPhlAn that cannot be run is reported rather than ignored"""

        with self.assertRaises(SystemExit):
            prescreen.verify_metaphlan_db_version(
                "vOct22_CHOCOPhlAnSGB_202403",
                os.path.join(self.tempdir, "does_not_exist"))


class TestMetaphlanCommandLine(unittest.TestCase):
    """
    Test the command line prescreen.alignment builds for MetaPhlAn
    """

    def setUp(self):
        logging.getLogger("humann.search.prescreen").addHandler(logging.NullHandler())

        self.tempdir = tempfile.mkdtemp(prefix="humann_mpa_cmd_")
        config.temp_dir = self.tempdir
        config.unnamed_temp_dir = self.tempdir
        config.file_basename = "test"
        config.threads = 1

        self.saved_opts = config.metaphlan_opts
        self.captured = []
        self.saved = (prescreen.verify_metaphlan_db_version,
                      prescreen.utilities.execute_command)
        prescreen.verify_metaphlan_db_version = lambda *a, **k: None
        prescreen.utilities.execute_command = \
            lambda exe, args, *a, **k: self.captured.append([str(arg) for arg in args])

    def tearDown(self):
        (prescreen.verify_metaphlan_db_version,
         prescreen.utilities.execute_command) = self.saved
        config.metaphlan_opts = self.saved_opts
        shutil.rmtree(self.tempdir, ignore_errors=True)

    def metaphlan_args(self):
        prescreen.alignment(cfg.demo_fastq)
        return self.captured[0]

    def test_database_index_is_pinned(self):
        """The database HUMAnN verified is the one MetaPhlAn is asked for"""

        args = self.metaphlan_args()
        self.assertIn("--index", args)
        self.assertIn(config.metaphlan_v4_db_index, args)

    def test_pinned_index_carries_the_mpa_prefix(self):
        """MetaPhlAn resolves --index straight to a filename, so mpa_ is required"""

        self.assertEqual(config.metaphlan_v4_db_index,
                         "mpa_" + config.metaphlan_v4_db_version)

    def test_user_supplied_index_is_not_overridden(self):
        """An index given through --metaphlan-options wins"""

        config.metaphlan_opts = ["-t", "rel_ab_w_read_stats", "-x", "mpa_custom"]
        args = self.metaphlan_args()
        self.assertEqual(args.count("-x") + args.count("--index"), 1)
        self.assertNotIn(config.metaphlan_v4_db_index, args)


if __name__ == "__main__":
    unittest.main()
