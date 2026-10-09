"""Focused coverage for the experimental MetaPhlAn version-check override."""
import contextlib
import io
import os
from pathlib import Path
import runpy
import sys
import tempfile
import unittest
from unittest import mock

from humann import config, humann
from humann.search import prescreen


class TestMetaphlanOverride(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        saved = dict(config.__dict__)
        self.addCleanup(config.__dict__.update, saved)
        self.input = os.path.join(self.temp.name, "reads.fastq")
        Path(self.input).write_text("@read\nACGT\n+\nIIII\n")
        self.output = os.path.join(self.temp.name, "output")

    def args(self, *options):
        with mock.patch.object(sys, "argv", ["humann", "-i", self.input, "-o", self.output] + list(options)):
            return humann.parse_arguments(sys.argv)

    def test_override_is_opt_in(self):
        self.assertFalse(self.args().override_metaphlan_check)
        self.assertTrue(self.args("--override-metaphlan-check").override_metaphlan_check)

    def test_override_logs_critical_and_sets_runtime_configuration(self):
        args = self.args("--override-metaphlan-check", "--log-level", "CRITICAL")
        stderr = io.StringIO()
        with mock.patch.object(humann.logging, "basicConfig"), \
             self.assertLogs("humann.humann", level="CRITICAL") as log, \
             contextlib.redirect_stdout(io.StringIO()), \
             contextlib.redirect_stderr(stderr):
            humann.update_configuration(args)
        self.assertTrue(config.override_metaphlan_check)
        self.assertIn("EXPERIMENTAL", log.output[0])
        self.assertIn("take no responsibility", log.output[0])
        self.assertIn("CRITICAL", stderr.getvalue())
        # Omission on a subsequent invocation must reset the opt-in state.
        with mock.patch.object(humann.logging, "basicConfig"), \
             mock.patch.object(humann.logger, "critical") as warning, \
             contextlib.redirect_stdout(io.StringIO()):
            humann.update_configuration(self.args())
        self.assertFalse(config.override_metaphlan_check)
        warning.assert_not_called()

    def test_prescreen_skips_only_version_probe_when_overridden(self):
        for override in (False, True):
            with mock.patch.object(config, "override_metaphlan_check", override), \
                 mock.patch.object(prescreen, "verify_metaphlan_db_version") as verify, \
                 mock.patch.object(prescreen.utilities, "return_exe_path", return_value="/fake/bin"), \
                 mock.patch.object(prescreen.utilities, "fasta_or_fastq", return_value="fastq"), \
                 mock.patch.object(prescreen.utilities, "name_temp_file", return_value="/fake/bowtie2.txt"), \
                 mock.patch.object(prescreen.utilities, "execute_command") as execute, \
                 contextlib.redirect_stdout(io.StringIO()):
                prescreen.alignment(self.input)
            if override:
                verify.assert_not_called()
            else:
                verify.assert_called_once_with(config.metaphlan_v4_db_version, "metaphlan")
            self.assertEqual(execute.call_args[0][0], "metaphlan")

    def test_default_prescreen_still_rejects_failed_version_check(self):
        with mock.patch.object(config, "override_metaphlan_check", False), \
             mock.patch.object(prescreen, "verify_metaphlan_db_version", side_effect=SystemExit("incompatible")), \
             mock.patch.object(prescreen.utilities, "execute_command") as execute:
            with self.assertRaises(SystemExit):
                prescreen.alignment(self.input)
        execute.assert_not_called()

    def requirement_settings(self):
        database = os.path.join(self.temp.name, "chocophlan")
        os.makedirs(database)
        Path(database, "SGB0001.ffn").write_text(">gene\nACGT\n")
        config.nucleotide_database = database
        config.bypass_nucleotide_index = False
        config.bypass_prescreen = False
        config.bypass_nucleotide_search = True
        config.bypass_translated_search = True
        config.output_format = "tsv"

    def test_precomputed_profile_does_not_require_metaphlan(self):
        self.requirement_settings()
        profile = os.path.join(self.temp.name, "profile.tsv")
        Path(profile).write_text("# precomputed profile\n")
        for options in ([], ["--override-metaphlan-check"]):
            args = self.args("--taxonomic-profile", profile, "--input-format", "fastq", *options)
            with mock.patch.object(humann.utilities, "file_exists_readable"), \
                 mock.patch.object(humann.utilities, "find_exe_in_path", return_value=False) as find:
                humann.check_requirements(args)
            find.assert_not_called()

    def test_running_prescreen_still_requires_executable_with_override(self):
        self.requirement_settings()
        args = self.args("--input-format", "fastq", "--override-metaphlan-check")
        with mock.patch.object(humann.utilities, "file_exists_readable"), \
             mock.patch.object(humann.utilities, "find_exe_in_path", return_value=False) as find:
            with self.assertRaisesRegex(SystemExit, "metaphlan executable"):
                humann.check_requirements(args)
        find.assert_called_once_with("metaphlan")

    def test_override_keeps_precomputed_profile_header_validation(self):
        config.override_metaphlan_check = True
        profile = os.path.join(self.temp.name, "invalid-profile.tsv")
        Path(profile).write_text("#" + config.metaphlan_v4_db_version + "\n# invalid columns\n")
        with mock.patch.object(prescreen.utilities, "name_temp_file", return_value="/unused/custom.ffn"), \
             self.assertLogs("humann.search.prescreen", level="ERROR"):
            with self.assertRaisesRegex(SystemExit, "expected columns"):
                prescreen.create_custom_database(self.temp.name, profile)


class TestInstallMetaphlanOverride(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import setuptools
        with mock.patch.object(setuptools, "setup"):
            state = runpy.run_path(str(Path(__file__).resolve().parents[2] / "setup.py"))
        cls.install_class = state["Install"]
        cls.install_globals = cls.install_class.run.__globals__

    def command(self, override):
        from setuptools import Distribution
        options = ["install"] + (["--override-metaphlan-check"] if override else [])
        dist = Distribution({"script_name": "setup.py", "script_args": options,
                             "cmdclass": {"install": self.install_class}})
        with contextlib.redirect_stdout(io.StringIO()):
            dist.parse_command_line()
        return dist.get_command_obj("install")

    def run_command(self, command, folder):
        command.install_lib = folder
        command.bypass_dependencies_install = True
        with mock.patch.dict(self.install_globals,
                             {"check_metaphlan_version": mock.Mock(), "install_minpath": mock.Mock()}), \
             mock.patch.object(self.install_globals["_install"], "run"), \
             contextlib.redirect_stdout(io.StringIO()):
            check = self.install_globals["check_metaphlan_version"]
            command.run()
            return check.call_count

    def test_default_install_checks_before_and_after_installation(self):
        command = self.command(False)
        self.assertFalse(command.override_metaphlan_check)
        with tempfile.TemporaryDirectory() as folder:
            self.assertEqual(self.run_command(command, folder), 2)

    def test_install_override_skips_both_checks_and_logs_critical(self):
        command = self.command(True)
        self.assertTrue(command.override_metaphlan_check)
        with tempfile.TemporaryDirectory() as folder, self.assertLogs(level="CRITICAL") as log:
            self.assertEqual(self.run_command(command, folder), 0)
        self.assertIn("EXPERIMENTAL", log.output[0])
        self.assertIn("take no responsibility", log.output[0])


if __name__ == "__main__":
    unittest.main()
