"""Nonfatal deprecation warnings for legacy pathway profiles on ingestion."""
import bz2
import contextlib
import gzip
import io
import sys
import tempfile
import unittest
from pathlib import Path
from types import ModuleType
from unittest import mock

from humann import accounting, config, store
from humann.tools import util, split_table, expand_cluster


class TestLegacyProfileWarning(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.lines = ['# Pathway\tsample\n', 'UNMAPPED\t0\n',
                      'UNINTEGRATED: old|g__A.s__a\t123.456789\n',
                      'UNINTEGRATED|g__B.s__b\t10\n', 'P1\t42\n']
        self.path = self.root / 'legacy.tsv'
        self.path.write_text(''.join(self.lines))

    def test_stream_warning_is_once_nonfatal_and_preserves_values(self):
        with self.assertLogs('humann.accounting', level='WARNING') as logged:
            result = list(accounting.warn_legacy_profile(iter(self.lines), str(self.path)))
        self.assertEqual(result, self.lines)
        self.assertEqual(len(logged.output), 1)
        warning = logged.output[0]
        for message in ('DEPRECATION', str(self.path), 'earlier version of HUMAnN',
                        'UNMAPPED', 'UNINTEGRATED', 'CPM-scale abundances',
                        'does not interrupt execution', 'legacy input values remain supported'):
            self.assertIn(message, warning)

    def test_either_legacy_name_triggers_warning_including_named_strata(self):
        for feature in ('UNMAPPED', 'UNINTEGRATED', 'UNMAPPED: old|taxon',
                        'UNINTEGRATED|taxon|nested'):
            with self.subTest(feature=feature), self.assertLogs('humann.accounting', level='WARNING') as logged:
                list(accounting.warn_legacy_profile([feature + '\t0\n'], 'input.tsv'))
            self.assertEqual(len(logged.output), 1)

    def test_modern_counts_and_biological_ids_do_not_trigger_warning(self):
        lines = ['# READS_INITIAL\t1000\n', '# READS_UNMAPPED\t200\n',
                 '# Pathway\tsample\n', 'READS_UNMAPPED\t200\n', 'UNGROUPED\t10\n',
                 'UNMAPPED_OTHER\t5\n', 'P1: UNINTEGRATED|UNMAPPED\t42\n']
        with mock.patch.object(accounting.logger, 'warning') as warning:
            self.assertEqual(list(accounting.warn_legacy_profile(lines, 'modern.tsv')), lines)
        warning.assert_not_called()

    def test_legacy_names_and_warning_text_inherit_config(self):
        names = ('OLD_UNMAPPED', 'OLD_UNINTEGRATED')
        with mock.patch.object(config, 'legacy_pathway_feature_names', names), \
             self.assertLogs('humann.accounting', level='WARNING') as logged:
            list(accounting.warn_legacy_profile(['OLD_UNINTEGRATED|taxon\t1\n'], 'custom.tsv'))
        self.assertIn('OLD_UNMAPPED and OLD_UNINTEGRATED', logged.output[0])
        with mock.patch.object(config, 'legacy_pathway_feature_names', names), \
             mock.patch.object(accounting.logger, 'warning') as warning:
            list(accounting.warn_legacy_profile(self.lines, 'custom.tsv'))
        warning.assert_not_called()

    def test_shared_table_reader_warns_once_per_ingestion_even_same_file(self):
        with self.assertLogs('humann.accounting', level='WARNING') as logged:
            for _ in range(2):
                table = util.Table(str(self.path))
                self.assertIn('UNMAPPED', table.rowheads)
                self.assertEqual(table.data[table.rowheads.index('UNINTEGRATED: old|g__A.s__a')], ['123.456789'])
        self.assertEqual(len(logged.output), 2)

    def test_shared_generator_reader_warns_for_non_table_utilities(self):
        with self.assertLogs('humann.accounting', level='WARNING') as logged:
            rows = list(util.process_gene_table_with_header(str(self.path)))
        self.assertEqual(len(logged.output), 1)
        self.assertEqual(rows[1], 'UNMAPPED\t0')

    def test_mapping_files_are_not_misclassified_as_profiles(self):
        with mock.patch.object(accounting.logger, 'warning') as warning:
            list(util.gzip_bzip2_biom_open_readlines(str(self.path)))
        warning.assert_not_called()

    def test_compressed_profiles_warn_once(self):
        for extension, opener in (('.gz', gzip.open), ('.bz2', bz2.open)):
            with self.subTest(extension=extension):
                path = str(self.root / ('legacy.tsv' + extension))
                with opener(path, 'wt') as handle:
                    handle.write(''.join(self.lines))
                with self.assertLogs('humann.accounting', level='WARNING') as logged:
                    table = util.Table(path)
                self.assertEqual(len(logged.output), 1)
                self.assertIn('UNMAPPED', table.rowheads)

    def test_stdin_profile_warns_once_and_keeps_rows(self):
        with mock.patch.object(sys, 'stdin', io.StringIO(''.join(self.lines))), \
             self.assertLogs('humann.accounting', level='WARNING') as logged:
            table = util.Table(None)
        self.assertEqual(len(logged.output), 1)
        self.assertIn('<STDIN>', logged.output[0])
        self.assertEqual(table.rowheads, ['UNMAPPED', 'UNINTEGRATED: old|g__A.s__a',
                                         'UNINTEGRATED|g__B.s__b', 'P1'])

    def test_core_gene_table_ingestion_warns_and_continues(self):
        scores = store.GeneScores()
        with self.assertLogs('humann.accounting', level='WARNING') as logged:
            scores.add_from_file(str(self.path))
        self.assertEqual(len(logged.output), 1)
        self.assertEqual(scores.scores_for_bug('all'), {'P1': 42})

    def test_direct_sample_splitter_warns_without_altering_legacy_values(self):
        output = self.root / 'split'; output.mkdir()
        with self.assertLogs('humann.accounting', level='WARNING') as logged:
            files = split_table.split_gene_table(str(self.path), str(output))
        self.assertEqual(len(logged.output), 1)
        self.assertIn('UNINTEGRATED: old|g__A.s__a\t123.456789', Path(files[0]).read_text())

    def test_cluster_selection_warns_and_carries_legacy_rows(self):
        output = str(self.root / 'cluster.tsv')
        with self.assertLogs('humann.accounting', level='WARNING') as logged:
            expand_cluster.write_genes([], {}, str(self.path), output, False)
        self.assertEqual(len(logged.output), 1)
        self.assertIn('UNINTEGRATED: old|g__A.s__a\t123.456789', Path(output).read_text())

    def test_biom_profile_warns_once_through_table_reader_or_conversion(self):
        # Fake only the BIOM loading API; exercise the real profile-reader paths.
        biom = ModuleType('biom')
        loaded = mock.Mock()
        loaded.ids.return_value = ['sample']
        loaded.metadata.return_value = {}
        loaded.to_tsv.return_value = '# Constructed from biom file\n' + ''.join(self.lines).rstrip('\n')
        biom.load_table = mock.Mock(return_value=loaded)
        path = self.root / 'legacy.biom'; path.write_bytes(b'fake-biom')
        with mock.patch.dict(sys.modules, {'biom': biom}):
            with self.assertLogs('humann.accounting', level='WARNING') as logged:
                table = util.Table(str(path))
            self.assertEqual(len(logged.output), 1)
            self.assertIn('UNMAPPED', table.rowheads)
            output = str(self.root / 'converted.tsv')
            with self.assertLogs('humann.accounting', level='WARNING') as logged:
                util.biom_to_tsv(str(path), output)
            self.assertEqual(len(logged.output), 1)
            self.assertIn('UNMAPPED\t0', Path(output).read_text())


if __name__ == '__main__':
    unittest.main()
