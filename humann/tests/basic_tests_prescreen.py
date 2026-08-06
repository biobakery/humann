"""
Tests for the taxonomic prescreen and custom ChocoPhlAn database creation
"""

import unittest
import contextlib
import io
import logging
import shutil
import tempfile

import cfg
import utils

from humann import config
from humann.search import prescreen


@contextlib.contextmanager
def captured_stdout():
    """Collect what the prescreen prints for the user"""

    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        yield buffer


# the DEMO ChocoPhlAn ships pangenomes for exactly these taxa
DEMO_CHOCOPHLAN_TAXA = ["SGB1871", "SGB2091", "SGB2301"]
# selected by the fixture profile, with no pangenome in any ChocoPhlAn build
# keyed on the vOct22 database
TAXA_WITHOUT_PANGENOMES = ["SGB99999", "EUK100861"]


class TestPrescreenProfileParsing(unittest.TestCase):
    """
    Test prescreen.create_custom_database against a MetaPhlAn 4.2 profile
    """

    def setUp(self):
        logging.getLogger("humann.search.prescreen").addHandler(logging.NullHandler())

        self.tempdir = tempfile.mkdtemp(prefix="humann_prescreen_test_")
        config.temp_dir = self.tempdir
        config.unnamed_temp_dir = self.tempdir
        config.file_basename = "test"
        config.bypass_prescreen = False
        config.prescreen_threshold = 0.5
        config.average_read_length = 150
        config.sgb_to_species_mapping = {}

    def tearDown(self):
        shutil.rmtree(self.tempdir, ignore_errors=True)

    def selected_taxa(self, profile):
        """Run the prescreen and return the taxa it selected"""

        prescreen.create_custom_database(cfg.demo_chocophlan_dir, profile)
        return set(config.sgb_to_species_mapping.keys())

    def test_vJan25_profile_is_accepted(self):
        """A MetaPhlAn 4.2 / vJan25 profile passes the header checks"""

        # create_custom_database exits on an unrecognised header
        self.selected_taxa(cfg.metaphlan_profile_vJan25)

    def test_unclassified_row_is_skipped(self):
        """The UNCLASSIFIED row is not parsed as a clade

        MetaPhlAn 4.2 emits UNCLASSIFIED by default, and its coverage field is a
        literal '-' that cannot be read as a number.
        """

        self.assertNotIn("UNCLASSIFIED", self.selected_taxa(cfg.metaphlan_profile_vJan25))

    def test_taxa_above_threshold_are_selected(self):
        """Taxa whose coverage clears the prescreen threshold are selected"""

        selected = self.selected_taxa(cfg.metaphlan_profile_vJan25)
        self.assertIn("SGB1871", selected)
        self.assertIn("SGB2091", selected)

    def test_taxa_below_threshold_are_excluded(self):
        """SGB2301 falls below prescreen_threshold and is not selected"""

        self.assertNotIn("SGB2301", self.selected_taxa(cfg.metaphlan_profile_vJan25))

    def test_species_name_is_tracked(self):
        """The s__ label is captured for naming stratified output"""

        self.selected_taxa(cfg.metaphlan_profile_vJan25)
        self.assertEqual(config.sgb_to_species_mapping["SGB1871"], "s__Bacteroides_ovatus")

    def test_higher_rank_rows_are_ignored(self):
        """Rows without a t__ level contribute nothing to the database"""

        for taxon in self.selected_taxa(cfg.metaphlan_profile_vJan25):
            self.assertNotIn("|", taxon)


class TestReportTaxaWithoutPangenomes(unittest.TestCase):
    """
    Test prescreen.report_taxa_without_pangenomes
    """

    def setUp(self):
        logging.getLogger("humann.search.prescreen").addHandler(logging.NullHandler())

    def test_silent_when_every_taxon_matched(self):
        """Nothing is reported when ChocoPhlAn covered the whole selection"""

        with captured_stdout() as output:
            prescreen.report_taxa_without_pangenomes(
                [], DEMO_CHOCOPHLAN_TAXA, {"SGB1871": 30.4})
        self.assertEqual(output.getvalue(), "")

    def test_missing_taxa_are_named(self):
        """Unmatched taxa appear in the warning"""

        with captured_stdout() as output:
            prescreen.report_taxa_without_pangenomes(
                TAXA_WITHOUT_PANGENOMES,
                DEMO_CHOCOPHLAN_TAXA + TAXA_WITHOUT_PANGENOMES,
                {"SGB99999": 9.8, "EUK100861": 0.7})
        for taxon in TAXA_WITHOUT_PANGENOMES:
            self.assertIn(taxon, output.getvalue())

    def test_missing_abundance_is_summed(self):
        """The warning states how much of the community was stranded"""

        with captured_stdout() as output:
            prescreen.report_taxa_without_pangenomes(
                TAXA_WITHOUT_PANGENOMES,
                DEMO_CHOCOPHLAN_TAXA + TAXA_WITHOUT_PANGENOMES,
                {"SGB99999": 9.8, "EUK100861": 0.7})
        self.assertIn("10.50%", output.getvalue())

    def test_long_lists_are_truncated(self):
        """A mismatched database strands thousands of taxa; the list is capped"""

        missing = ["SGB" + str(number) for number in range(500)]
        with captured_stdout() as output:
            prescreen.report_taxa_without_pangenomes(missing, missing, {})
        printed = output.getvalue()
        self.assertIn("500 of the 500", printed)
        self.assertIn("and " + str(500 - prescreen.MAX_REPORTED_MISSING_TAXA) + " more",
                      printed)


class TestCustomDatabaseReporting(unittest.TestCase):
    """
    Test that create_custom_database reports its unmatched selections
    """

    def setUp(self):
        logging.getLogger("humann.search.prescreen").addHandler(logging.NullHandler())

        self.tempdir = tempfile.mkdtemp(prefix="humann_prescreen_report_")
        config.temp_dir = self.tempdir
        config.unnamed_temp_dir = self.tempdir
        config.file_basename = "test"
        config.bypass_prescreen = False
        config.prescreen_threshold = 0.5
        config.average_read_length = 150
        config.sgb_to_species_mapping = {}

    def tearDown(self):
        shutil.rmtree(self.tempdir, ignore_errors=True)

    def test_taxa_without_pangenomes_are_reported(self):
        """Selected taxa the DEMO ChocoPhlAn cannot serve are named on stdout

        The fixture selects four taxa; the DEMO database holds pangenomes for
        two of them. Before this was reported the run printed "Total species
        selected from prescreen: 4" and quietly built a database from two.
        """

        with captured_stdout() as output:
            prescreen.create_custom_database(
                cfg.demo_chocophlan_dir, cfg.metaphlan_profile_vJan25)

        printed = output.getvalue()
        self.assertIn("no pangenome in the ChocoPhlAn database", printed)
        for taxon in TAXA_WITHOUT_PANGENOMES:
            self.assertIn(taxon, printed)

    def test_matched_taxa_are_not_reported_as_missing(self):
        """Taxa the DEMO ChocoPhlAn does serve stay out of the warning"""

        with captured_stdout() as output:
            prescreen.create_custom_database(
                cfg.demo_chocophlan_dir, cfg.metaphlan_profile_vJan25)

        warning = [line for line in output.getvalue().splitlines()
                   if line.startswith("WARNING:")]
        self.assertEqual(len(warning), 1)
        self.assertNotIn("SGB1871", warning[0])
        self.assertNotIn("SGB2091", warning[0])


if __name__ == "__main__":
    unittest.main()
