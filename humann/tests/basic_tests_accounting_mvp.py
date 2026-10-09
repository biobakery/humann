"""Regression coverage for mapped-only CPM and protected raw-count metadata."""
import contextlib
import importlib
import importlib.util
import io
import sys
import runpy
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace, ModuleType
from unittest import mock

from humann import accounting, config, store
from humann.quantify import families, modules
from humann.tools import util, renorm_table, rename_table, join_tables, split_table
from humann.tools import split_stratified_table, genefamilies_genus_level, expand_cluster
from humann.tools import merge_abundance, rna_dna_norm, reduce_table

# Installed mapping databases are not needed for these custom-mapping tests.
with mock.patch('os.listdir', return_value=['map_ko_uniref90.txt.gz']):
    regroup_table = importlib.import_module('humann.tools.regroup_table')
with mock.patch.object(sys, 'argv', ['infer_taxonomy', '--devdb', 'unused']):
    infer_taxonomy = importlib.import_module('humann.tools.infer_taxonomy')


class TestAccountingMVP(unittest.TestCase):
    def setUp(self):
        saved = dict(config.__dict__)
        self.addCleanup(config.__dict__.update, saved)
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        config.genefamilies_file = str(self.root / 'genes.tsv')
        config.reactions_file = str(self.root / 'reactions.tsv')
        config.pathabundance_file = str(self.root / 'pathways.tsv')
        config.pathcoverage_file = str(self.root / 'coverage.tsv')
        config.file_basename = 'sample'
        config.output_format = 'tsv'
        config.count_normalization = 'Adjusted CPMs'
        config.remove_stratified_output = False
        config.minpath_toggle = 'off'
        self.counts = {config.initial_reads_name: ['10000000'], config.unmapped_gene_name: ['2000000']}

    def table(self, rows, name='input.tsv', counts=True):
        path = self.root / name
        comments = list(accounting.comment_lines(self.counts)) if counts else []
        path.write_text('\n'.join(comments + ['# Feature\tsample'] + rows) + '\n')
        return util.Table(str(path))

    def scores(self, unmapped=2000000):
        alignments = store.Alignments()
        alignments.add('G1', 1000, 'read1', 100, 'g__A.s__a', read_length=101)
        alignments.add('G2', 1900, 'read2', 100, 'g__B.s__b', read_length=101)
        scores = store.GeneScores()
        families.gene_families(alignments, scores, unmapped, initial_read_count=10000000)
        return scores

    def test_mapped_cpm_denominator_and_downstream_scores(self):
        scores = self.scores()
        self.assertAlmostEqual(scores.get_score('all', 'G1'), 2e6 / 3)
        self.assertAlmostEqual(scores.get_score('all', 'G2'), 1e6 / 3)
        self.assertAlmostEqual(sum(scores.scores_for_bug('all').values()), 1e6)
        self.assertEqual(scores.get_score('g__A.s__a', 'G1'), scores.get_score('all', 'G1'))
        table = util.Table(config.genefamilies_file)
        self.assertEqual(table.read_counts, self.counts)
        self.assertFalse(any(util.is_special(h) for h in table.rowheads))
        self.assertAlmostEqual(float(table.data[table.rowheads.index('G1')][0]), 2e6 / 3, places=6)
        other = self.scores(unmapped=9000000)
        self.assertEqual(other.scores_for_bug('all'), scores.scores_for_bug('all'))

    def test_empty_mapped_fraction_retains_counts_without_fake_features(self):
        families.gene_families(store.Alignments(), store.GeneScores(), 10, initial_read_count=10)
        table = util.Table(config.genefamilies_file)
        self.assertEqual(table.rowheads, [])
        self.assertEqual(table.read_counts[config.unmapped_gene_name], ['10'])

    def test_rpk_mode_retains_existing_length_convention(self):
        config.count_normalization = 'RPKs'
        scores = self.scores()
        self.assertEqual(scores.get_score('all', 'G1'), 1)
        self.assertAlmostEqual(scores.get_score('all', 'G2'), 1 / 1.9)

    def test_counts_logged_independently_of_header_rendering(self):
        with self.assertLogs('humann.accounting', level='INFO') as logged:
            accounting.log_counts({config.initial_reads_name: 10, config.unmapped_gene_name: 2})
        self.assertIn('READS_INITIAL\t10', logged.output[0])
        self.assertIn('READS_UNMAPPED\t2', logged.output[1])

    def test_reserved_names_are_config_driven(self):
        with mock.patch.object(config, 'special_feature_names', ('CUSTOM_RESERVED',)):
            self.assertTrue(util.is_special('CUSTOM_RESERVED: label|taxon|nested'))
            self.assertFalse(util.is_special('UNMAPPED'))

    def test_old_pathway_accounting_never_becomes_gene_assignments(self):
        table = self.table(['UNMAPPED\t900', 'UNINTEGRATED|g__A.s__a\t50', 'UNGROUPED\t5', 'G1\t100'])
        scores = store.GeneScores()
        self.assertEqual(scores.add_from_file(str(self.root / 'input.tsv')), 2000000)
        self.assertEqual(scores.scores_for_bug('all'), {'G1': 100})
        self.assertEqual(scores.read_counts[config.initial_reads_name], 10000000)

    def test_legacy_unmapped_read_count_fallback_and_unknown_initial(self):
        self.table(['READS_UNMAPPED\t123', 'G1\t100'], counts=False)
        scores = store.GeneScores()
        self.assertEqual(scores.add_from_file(str(self.root / 'input.tsv')), 123)
        self.assertIsNone(scores.read_counts[config.initial_reads_name])

    def test_legacy_pathway_unmapped_does_not_fabricate_a_read_count(self):
        self.table(['UNMAPPED\t123', 'G1\t100'], counts=False)
        scores = store.GeneScores()
        self.assertIsNone(scores.add_from_file(str(self.root / 'input.tsv')))

    def test_regroup_preserves_accounting_and_allows_overlapping_groups(self):
        table = self.table(['UNMAPPED: original\t123.456789', 'UNINTEGRATED|g__A.s__a\t0', 'G1\t600000', 'G2\t400000'])
        regroup_table.regroup(table, {'UNMAPPED': ['BAD'], 'G1': ['K1', 'K2']}, 'mean', 2, ungrouped=True)
        values = dict(zip(table.rowheads, table.data))
        self.assertNotIn('BAD', values)
        self.assertEqual(values['UNMAPPED: original'], ['123.456789'])
        self.assertEqual(values['UNINTEGRATED|g__A.s__a'], ['0'])
        self.assertEqual(values['K1'], [600000])
        self.assertEqual(values['K2'], [600000])
        self.assertEqual(values[config.ungrouped_reaction_name], [400000])
        self.assertEqual(table.read_counts, self.counts)

    def test_regroup_residual_is_summed_even_under_mean(self):
        table = self.table(['G1\t10', 'G2\t20'])
        regroup_table.regroup(table, {}, 'mean', 2, ungrouped=True)
        self.assertEqual(table.data, [[30]])
        with self.assertRaises(ValueError):
            regroup_table.regroup(self.table(['G1\t10']), {'G1': ['READS_UNMAPPED']}, 'sum', None)

    def test_reactions_derive_from_cpm_and_keep_count_comments(self):
        scores = self.scores()
        reactions = mock.Mock()
        reactions.reaction_list.return_value = ['R1', 'R2']
        reactions.find_genes.return_value = ['G1']
        pathways = mock.Mock()
        pathways.find_pathways.return_value = ['P1']
        modules.identify_reactions_and_pathways(scores, reactions, pathways, 2000000)
        table = util.Table(config.reactions_file)
        values = dict(zip(table.rowheads, table.data))
        self.assertAlmostEqual(float(values['R1'][0]), 2e6 / 3, places=6)
        self.assertAlmostEqual(float(values['R2'][0]), 2e6 / 3, places=6)
        self.assertAlmostEqual(float(values[config.ungrouped_reaction_name][0]), 1e6 / 3, places=6)
        self.assertNotIn('UNMAPPED', values)
        self.assertEqual(table.read_counts, self.counts)

    def test_no_reaction_mapping_does_not_classify_assigned_genes_as_ungrouped(self):
        scores = self.scores()
        pathways = mock.Mock()
        pathways.find_pathways.return_value = ['P1']
        modules.identify_reactions_and_pathways(scores, None, pathways, 2000000)
        values = dict(zip(util.Table(config.reactions_file).rowheads, util.Table(config.reactions_file).data))
        self.assertEqual(float(values[config.ungrouped_reaction_name][0]), 0)
        self.assertIn('G1', values)

    def test_pathway_output_has_no_legacy_residuals(self):
        scores = self.scores()
        abundance = store.Pathways()
        abundance.add('all', 'P1', 42)
        database = mock.Mock()
        database.get_pathways_names.return_value = {}
        with mock.patch.object(modules, 'compute_pathways_abundance', return_value=(abundance, {})):
            modules.compute_pathways_abundance_and_coverage(scores, None, store.PathwaysAndReactions(), database, 2000000)
        table = util.Table(config.pathabundance_file)
        self.assertEqual(table.rowheads, ['P1'])
        self.assertEqual(float(table.data[0][0]), 42)
        self.assertEqual(table.read_counts, self.counts)

    def test_real_pathway_computation_uses_normalized_gene_scores(self):
        scores = self.scores()
        reactions = store.ReactionsDatabase()
        reactions.add_reactions({'R1': ['G1'], 'R2': ['G2']})
        pathways = store.PathwaysDatabase()
        pathways.add_pathway_structure('P1', ' R1 R2 ')
        config.gap_fill_toggle = 'off'
        linked = modules.identify_reactions_and_pathways(scores, reactions, pathways, 2000000)
        modules.compute_pathways_abundance_and_coverage(scores, reactions, linked, pathways, 2000000)
        table = util.Table(config.pathabundance_file)
        values = dict(zip(table.rowheads, table.data))
        # Harmonic mean of 2/3 million and 1/3 million, with no pathway re-TSS.
        self.assertAlmostEqual(float(values['P1'][0]), 4e6 / 9, places=6)
        self.assertFalse(any(util.is_special(h) for h in table.rowheads))
        self.assertEqual(table.read_counts, self.counts)

    @unittest.skipUnless(all(importlib.util.find_spec(name) for name in ('biom', 'h5py', 'numpy')),
                         'BIOM round-trip dependencies are not installed')
    def test_biom_round_trip_keeps_counts_as_sample_metadata(self):
        table = self.table(['UNMAPPED\t123.456789', 'G1\t10'])
        biom_path = str(self.root / 'counts.biom')
        util.tsv_to_biom(str(self.root / 'input.tsv'), biom_path)
        restored = util.Table(biom_path)
        self.assertEqual(restored.read_counts, self.counts)
        self.assertNotIn(config.initial_reads_name, restored.rowheads)
        tsv_path = str(self.root / 'restored.tsv')
        util.biom_to_tsv(biom_path, tsv_path)
        self.assertEqual(util.Table(tsv_path).read_counts, self.counts)

    def test_renorm_preserves_legacy_rows_outside_denominator(self):
        table = self.table(['UNMAPPED\t999', 'UNINTEGRATED: old|g__A.s__a\t123.456789', 'G1\t2', 'G2\t3'])
        renorm_table.normalize(table)
        values = dict(zip(table.rowheads, table.data))
        self.assertEqual(values['UNMAPPED'], ['999'])
        self.assertEqual(values['UNINTEGRATED: old|g__A.s__a'], ['123.456789'])
        self.assertAlmostEqual(float(values['G1'][0]), 400000)
        self.assertAlmostEqual(float(values['G2'][0]), 600000)
        self.assertEqual(table.read_counts, self.counts)

    def test_join_and_both_split_modes_preserve_counts(self):
        self.table(['G1\t10', 'G1|g__A.s__a\t10'], name='one.tsv')
        self.table(['G2\t20'], name='two.tsv', counts=False)
        joined = str(self.root / 'joined.tsv')
        join_tables.join_gene_tables([str(self.root / 'one.tsv'), str(self.root / 'two.tsv')], joined)
        self.assertEqual(util.Table(joined).read_counts[config.initial_reads_name], ['10000000', 'NA'])
        out = self.root / 'split'; out.mkdir()
        files = split_table.split_gene_table(joined, str(out))
        self.assertEqual(util.Table(files[0]).read_counts, self.counts)
        self.assertEqual(util.Table(files[1]).read_counts[config.unmapped_gene_name], ['NA'])
        strata_files = split_stratified_table.split_table(str(self.root / 'one.tsv'), str(out))
        for path in strata_files:
            self.assertEqual(util.Table(path).read_counts, self.counts)

    def test_rename_genus_and_cluster_preserve_accounting(self):
        source = self.table(['UNINTEGRATED: old|g__A.s__a\t123.456789', 'G1\t10', 'G1|g__A.s__a\t10'])
        rename_table.rename(source, {'UNINTEGRATED': {'wrong': 1}, 'G1': {'name': 1}})
        self.assertIn('UNINTEGRATED: old|g__A.s__a', source.rowheads)
        out = str(self.root / 'genus.tsv')
        genefamilies_genus_level.create_table(str(self.root / 'input.tsv'), out)
        self.assertIn('UNINTEGRATED: old|g__A.s__a', util.Table(out).rowheads)
        cluster = str(self.root / 'cluster.tsv')
        expand_cluster.write_genes(['U50'], {'G1': ['U50']}, str(self.root / 'input.tsv'), cluster, False)
        self.assertIn('UNINTEGRATED: old|g__A.s__a', util.Table(cluster).rowheads)
        self.assertEqual(util.Table(cluster).read_counts, self.counts)

    def test_compressed_table_and_cluster_keep_counts(self):
        table = self.table(['UNMAPPED|taxon|nested\t123.456789', 'G1\t10'])
        for suffix in ('.gz', '.bz2'):
            path = str(self.root / ('compressed.tsv' + suffix))
            table.write(path)
            restored = util.Table(path)
            self.assertEqual(restored.read_counts, self.counts)
            self.assertEqual(restored.data, table.data)
        cluster = str(self.root / 'cluster.tsv')
        expand_cluster.write_genes(['U50'], {'G1': ['U50']},
                                   str(self.root / 'compressed.tsv.gz'), cluster, False)
        self.assertEqual(util.Table(cluster).read_counts, self.counts)
        regroup_table.regroup(table, {'G1': ['K1']}, 'sum', 2)
        self.assertIn('UNMAPPED|taxon|nested', table.rowheads)

    def test_merge_does_not_unpack_historical_accounting_rows(self):
        self.table(['UNMAPPED\t0', 'UNINTEGRATED|g__A.s__a\t123.456789', 'P1\t10'])
        out = str(self.root / 'merged.tsv')
        merge_abundance.merge_abundances({}, {'UNMAPPED': ['G1']}, str(self.root / 'input.tsv'), out, {}, False)
        table = util.Table(out)
        self.assertIn('UNMAPPED', table.rowheads)
        self.assertIn('UNINTEGRATED|g__A.s__a', table.rowheads)
        self.assertEqual(table.read_counts, self.counts)

    def test_taxonomy_inference_cannot_reassign_protected_rows(self):
        self.table(['UNINTEGRATED: old|g__A.s__a\t123.456789', 'G1\t10'])
        out = str(self.root / 'tax.tsv')
        args = SimpleNamespace(input=str(self.root / 'input.tsv'), output=out, devdb='unused', database=None,
                               level='Genus', lca_choice='humann_lca', threshold=0, mode=infer_taxonomy.c_tmode)
        with mock.patch.object(infer_taxonomy, 'get_args', return_value=args), \
             mock.patch.object(infer_taxonomy, 'build_taxmap', return_value={'G1': 'g__B', 'UNINTEGRATED': 'g__Wrong'}):
            infer_taxonomy.main()
        table = util.Table(out)
        self.assertIn('UNINTEGRATED: old|g__A.s__a', table.rowheads)
        self.assertNotIn('UNINTEGRATED|g__Wrong', table.rowheads)
        self.assertEqual(table.read_counts, self.counts)

    def test_table_reduction_does_not_invent_a_read_library(self):
        self.table(['UNMAPPED\t123.456789', 'G1\t10'])
        out = str(self.root / 'reduced.tsv')
        reduce_table.reduce_table('max', str(self.root / 'input.tsv'), out, False, None)
        self.assertEqual(util.Table(out).read_counts, self.counts)
        joined = self.root / 'joined.tsv'
        joined.write_text('# READS_INITIAL\t10\t20\n# READS_UNMAPPED\t1\t2\n# Feature\tA\tB\nG1\t10\t20\n')
        reduce_table.reduce_table('max', str(joined), out, False, None)
        self.assertEqual(util.Table(out).read_counts,
                         {name: ['NA'] for name in config.read_count_names})

    @unittest.skipUnless(all(importlib.util.find_spec(name) for name in ('matplotlib', 'scipy', 'numpy')),
                         'Plotting dependencies are not installed')
    def test_barplot_reads_the_sample_header_after_count_comments(self):
        barplot = importlib.import_module('humann.tools.humann_barplot')
        self.table(['G1\t10', 'G1|g__A.s__a\t10'])
        table = barplot.BarplotTable(str(self.root / 'input.tsv'), focal_feature='G1')
        self.assertEqual(table.colheads, ['sample'])
        self.assertEqual(table.data.tolist(), [[10]])

    def test_old_taxonomy_expander_preserves_protected_strata(self):
        self.table(['UNINTEGRATED|g__A.s__a\t123.456789', 'G1\t10'])
        out = str(self.root / 'expanded.tsv')
        ete = ModuleType('ete2')
        taxonomy = mock.Mock()
        ete.NCBITaxa = mock.Mock(return_value=taxonomy)
        argv = ['humann_expand_results_taxonomy', '--input', str(self.root / 'input.tsv'), '--output', out]
        with mock.patch.dict(sys.modules, {'ete2': ete}), mock.patch.object(sys, 'argv', argv), \
             contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(SystemExit) as result:
                runpy.run_module('humann.tools.humann_expand_results_taxonomy', run_name='__main__')
        self.assertEqual(result.exception.code, 0)
        taxonomy.get_name_translator.assert_not_called()
        self.assertEqual(Path(out).read_text(), (self.root / 'input.tsv').read_text())

    def test_rna_dna_accounting_bypasses_smoothing_and_expression(self):
        self.table(['UNMAPPED\t123.456789', 'G1\t10', 'G2\t20'], name='dna.tsv')
        self.table(['UNMAPPED\t987.654321', 'G1\t20', 'G2\t40'], name='rna.tsv')
        for method in ('laplace', 'witten_bell'):
            args = SimpleNamespace(input_dna=str(self.root / 'dna.tsv'), input_rna=str(self.root / 'rna.tsv'),
                                   output_basename=str(self.root / method), method=method, log_transform=False, log_base=2)
            with mock.patch.object(rna_dna_norm, 'get_args', return_value=args):
                rna_dna_norm.main()
            for extension, expected in ((rna_dna_norm.c_new_dna_extension, '123.456789'),
                                        (rna_dna_norm.c_new_rna_extension, '987.654321'),
                                        (rna_dna_norm.c_norm_rna_extension, '987.654321')):
                table = util.Table(args.output_basename + extension)
                self.assertEqual(table.data[table.rowheads.index('UNMAPPED')], [expected])
                self.assertEqual(table.read_counts, self.counts)
            expression = util.Table(args.output_basename + rna_dna_norm.c_norm_rna_extension)
            self.assertAlmostEqual(float(expression.data[expression.rowheads.index('G1')][0]),
                                   1 if method == 'laplace' else 32 / 31, places=5)


if __name__ == '__main__':
    unittest.main()
