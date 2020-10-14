import unittest
import tempfile
import os

import cfg
import utils

class TestFunctionalHumannTools(unittest.TestCase):
    """
    Test humann.tools
    """

    def test_humann_join_tables_tsv(self):
        """
        Test joining tsv files with humann_join_tables
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann_temp")
        
        # join the files
        utils.run_command(["humann_join_tables","--input",
                           cfg.data_folder,"--output",new_file,"--file_name",
                           cfg.multi_sample_genefamilies_split_basename,"--verbose"])
        
        # check the joined file is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.multi_sample_genefamilies))

        # remove the temp file
        utils.remove_temp_file(new_file)
 
    def test_humann_split_tables_tsv(self):
        """
        Test splitting a tsv file with humann_split_tables
        """

        input_file=cfg.multi_sample_genefamilies

        # create a temp directory
        temp_directory=utils.create_temp_folder("split_tables_tsv")

        # split the file
        utils.run_command(["humann_split_table","--input", input_file,
                           "--output",temp_directory,"--verbose"])

        # test the split files are as expected
        output_files=os.listdir(temp_directory)

        # sort the output files
        file_pairs=[]
        for file in output_files:
            filebasename=os.path.basename(file)
            # get the sample number for the file
            file=os.path.join(temp_directory,file)
            if filebasename[-1] == 1:
                file_pairs.append([file,cfg.multi_sample_genefamilies_split1])
            elif filebasename[-1] == 2:
                file_pairs.append([file,cfg.multi_sample_genefamilies_split2])

        for temp_file, file in file_pairs:
            self.assertTrue(utils.files_almost_equal(temp_file, file))

        # remove the temp folder
        utils.remove_temp_folder(temp_directory)
        

    def test_humann_regroup_table_uniref50_rxn_tsv(self):
        """
        Test regrouping the tsv file with humann_regroup_table
        Test with uniref50 to reactions mappings
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann_temp")
        
        # run the command
        utils.run_command(["humann_regroup_table","--input",cfg.regroup_input,"--output",
                           new_file,"--groups","uniref50_rxn"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.regroup_rxn_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann_regroup_table_uniref50_rxn_tsv_mean(self):
        """
        Test regrouping the tsv file with humann_regroup_table
        Test with uniref50 to reactions mappings
        Test with the mean instead of sum output
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann_temp")
        
        # run the command
        utils.run_command(["humann_regroup_table","--input",cfg.regroup_input,"--output",
                           new_file,"--groups","uniref50_rxn","--function","mean"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.regroup_rxn_mean_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann_regroup_table_custom_grouping_tsv(self):
        """
        Test regrouping the tsv file with humann_regroup_table
        Test with custom mappings
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann_temp")
        
        # run the command
        utils.run_command(["humann_regroup_table","--input",cfg.regroup_custom_input,"--output",
                           new_file,"--custom",cfg.regroup_custom_groups])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.regroup_custom_groups_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
       
    def test_humann_rename_table_uniref50_tsv(self):
        """
        Test renaming the tsv file entries with humann_rename_table
        Test with uniref50 names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann_temp")
        
        # run the command
        utils.run_command(["humann_rename_table","--input",cfg.rename_input,"--output",
                           new_file,"--names","uniref50"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.rename_uniref50_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann_rename_table_ko_tsv(self):
        """
        Test renaming the tsv file entries with humann_rename_table
        Test with ko names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann_temp")
        
        # run the command
        utils.run_command(["humann_rename_table","--input",cfg.rename_ko_input,"--output",
                           new_file,"--names","kegg-orthology"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.rename_ko_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann_rename_table_ec_tsv(self):
        """
        Test renaming the tsv file entries with humann_rename_table
        Test with ec names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann_temp")
        
        # run the command
        utils.run_command(["humann_rename_table","--input",cfg.rename_ec_input,"--output",
                           new_file,"--names","ec"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.rename_ec_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann_rename_table_rxn_tsv(self):
        """
        Test renaming the tsv file entries with humann_rename_table
        Test with rxn names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann_temp")
        
        # run the command
        utils.run_command(["humann_rename_table","--input",cfg.rename_rxn_input,"--output",
                           new_file,"--names","metacyc-rxn"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.rename_rxn_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann_rename_table_pathways_tsv(self):
        """
        Test renaming the tsv file entries with humann_rename_table
        Test with pathways names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann_temp")
        
        # run the command
        utils.run_command(["humann_rename_table","--input",cfg.rename_pathway_input,"--output",
                           new_file,"--names","metacyc-pwy"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.rename_pathway_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann_rename_table_custom_tsv(self):
        """
        Test renaming the tsv file entries with humann_rename_table
        Test with custom names file
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann_temp")
        
        # run the command
        utils.run_command(["humann_rename_table","--input",cfg.rename_input,"--output",
                           new_file,"--custom",cfg.rename_custom_mapping])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.rename_custom_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
            
            
    def test_humann_renorm_table_cpm_tsv(self):
        """
        Test renorm the tsv file entries with humann_renorm_table
        Test with cpm
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann_temp")
        
        # run the command
        utils.run_command(["humann_renorm_table","--input",cfg.renorm_input,"--output",
                           new_file,"--units","cpm"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.renorm_cpm_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann_renorm_table_relab_tsv(self):
        """
        Test renorm the tsv file entries with humann_renorm_table
        Test with relab
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann_temp")
        
        # run the command
        utils.run_command(["humann_renorm_table","--input",cfg.renorm_input,"--output",
                           new_file,"--units","relab"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.renorm_relab_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann_unpack_pathways_tsv(self):
        """
        Test the tsv gene families and pathway abundance file entries with humann_unpack_pathways
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann_temp")
        
        # run the command
        utils.run_command(["humann_unpack_pathways","--input-genes",cfg.merge_abundance_genefamilies_input,
                           "--input-pathways",cfg.merge_abundance_pathways_input,"--output",
                           new_file])
        
        # check the output file is as expected
        # allow for varying precision in the calculations with almost equal
        self.assertTrue(utils.files_almost_equal(new_file, cfg.merge_abundance_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann_unpack_pathways_remove_taxonomy_tsv(self):
        """
        Test the tsv gene families and pathway abundance file entries with humann_unpack_pathways
        Test with the remove taxonomy option which stratifies by pathway then gene instead of
        stratifying by pathway, taxonomy, then gene
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann_temp")
        
        # run the command
        utils.run_command(["humann_unpack_pathways","--input-genes",cfg.merge_abundance_genefamilies_input,
                           "--input-pathways",cfg.merge_abundance_pathways_input,"--output",
                           new_file,"--remove-taxonomy"])
        
        # check the output file is as expected
        # allow for varying precision in the calculations with almost equal
        self.assertTrue(utils.files_almost_equal(new_file, cfg.merge_abundance_remove_taxonomy_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
