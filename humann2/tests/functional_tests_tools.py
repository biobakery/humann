import unittest
import tempfile
import os

import cfg
import utils

class TestFunctionalHumann2Tools(unittest.TestCase):
    """
    Test humann2.tools
    """

    def test_humann2_join_tables_tsv(self):
        """
        Test joining tsv files with humann2_join_tables
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # join the files
        utils.run_command(["humann2_join_tables","--input",
                           cfg.data_folder,"--output",new_file,"--file_name",
                           cfg.multi_sample_genefamilies_split_basename,"--verbose"])
        
        # check the joined file is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.multi_sample_genefamilies))

        # remove the temp file
        utils.remove_temp_file(new_file)
 
    def test_humann2_split_tables_tsv(self):
        """
        Test splitting a tsv file with humann2_split_tables
        """

        input_file=cfg.multi_sample_genefamilies

        # create a temp directory
        temp_directory=utils.create_temp_folder("split_tables_tsv")

        # split the file
        utils.run_command(["humann2_split_table","--input", input_file,
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
        

    def test_humann2_regroup_table_uniref50_rxn_tsv(self):
        """
        Test regrouping the tsv file with humann2_regroup_table
        Test with uniref50 to reactions mappings
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_regroup_table","--input",cfg.regroup_input,"--output",
                           new_file,"--groups","uniref50_rxn"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.regroup_rxn_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_regroup_table_uniref50_rxn_tsv_mean(self):
        """
        Test regrouping the tsv file with humann2_regroup_table
        Test with uniref50 to reactions mappings
        Test with the mean instead of sum output
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_regroup_table","--input",cfg.regroup_input,"--output",
                           new_file,"--groups","uniref50_rxn","--function","mean"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.regroup_rxn_mean_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_regroup_table_custom_grouping_tsv(self):
        """
        Test regrouping the tsv file with humann2_regroup_table
        Test with custom mappings
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_regroup_table","--input",cfg.regroup_custom_input,"--output",
                           new_file,"--custom",cfg.regroup_custom_groups])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.regroup_custom_groups_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
       
    def test_humann2_rename_table_uniref50_tsv(self):
        """
        Test renaming the tsv file entries with humann2_rename_table
        Test with uniref50 names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_rename_table","--input",cfg.rename_input,"--output",
                           new_file,"--names","uniref50"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.rename_uniref50_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_rename_table_ko_tsv(self):
        """
        Test renaming the tsv file entries with humann2_rename_table
        Test with ko names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_rename_table","--input",cfg.rename_ko_input,"--output",
                           new_file,"--names","kegg-orthology"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.rename_ko_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_rename_table_ec_tsv(self):
        """
        Test renaming the tsv file entries with humann2_rename_table
        Test with ec names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_rename_table","--input",cfg.rename_ec_input,"--output",
                           new_file,"--names","ec"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.rename_ec_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_rename_table_rxn_tsv(self):
        """
        Test renaming the tsv file entries with humann2_rename_table
        Test with rxn names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_rename_table","--input",cfg.rename_rxn_input,"--output",
                           new_file,"--names","metacyc-rxn"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.rename_rxn_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_rename_table_pathways_tsv(self):
        """
        Test renaming the tsv file entries with humann2_rename_table
        Test with pathways names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_rename_table","--input",cfg.rename_pathway_input,"--output",
                           new_file,"--names","metacyc-pwy"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.rename_pathway_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_rename_table_custom_tsv(self):
        """
        Test renaming the tsv file entries with humann2_rename_table
        Test with custom names file
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_rename_table","--input",cfg.rename_input,"--output",
                           new_file,"--custom",cfg.rename_custom_mapping])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.rename_custom_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
            
            
    def test_humann2_renorm_table_cpm_tsv(self):
        """
        Test renorm the tsv file entries with humann2_renorm_table
        Test with cpm
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_renorm_table","--input",cfg.renorm_input,"--output",
                           new_file,"--units","cpm"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.renorm_cpm_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_renorm_table_relab_tsv(self):
        """
        Test renorm the tsv file entries with humann2_renorm_table
        Test with relab
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_renorm_table","--input",cfg.renorm_input,"--output",
                           new_file,"--units","relab"])
        
        # check the output is as expected
        self.assertTrue(utils.files_almost_equal(new_file, cfg.renorm_relab_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_rna_dna_norm_laplace_tsv(self):
        """
        Test norm the tsv file entries from dna and rna input files with humann2_rna_dna_norm_table
        Test with laplace
        """
        
        # create a temp folder
        tempdir=utils.create_temp_folder("rna_dna_norm_laplace")
        output_basename=os.path.join(tempdir,"rna_dna_norm")
        
        # run the command
        utils.run_command(["humann2_rna_dna_norm","--input_dna",cfg.rna_dna_norm_dna_input,
                           "--input_rna",cfg.rna_dna_norm_rna_input,"--output_basename",
                           output_basename,"--method","laplace"])
        
        # check the output files are as expected
        for file_extension, expected_output_file in zip(cfg.rna_dna_norm_file_names, cfg.rna_dna_norm_laplace_output_files):
            self.assertTrue(utils.files_almost_equal(output_basename+file_extension, expected_output_file))

        # remove the temp file
        utils.remove_temp_folder(tempdir)
        
    def test_humann2_rna_dna_norm_witten_bell_tsv(self):
        """
        Test norm the tsv file entries from dna and rna input files with humann2_rna_dna_norm_table
        Test with witten bell
        """
        
        # create a temp folder
        tempdir=utils.create_temp_folder("rna_dna_norm_witten_bell")
        output_basename=os.path.join(tempdir,"rna_dna_norm")
        
        # run the command
        utils.run_command(["humann2_rna_dna_norm","--input_dna",cfg.rna_dna_norm_dna_input,
                           "--input_rna",cfg.rna_dna_norm_rna_input,"--output_basename",
                           output_basename,"--method","witten_bell"])
        
        # check the output files are as expected
        for file_extension, expected_output_file in zip(cfg.rna_dna_norm_file_names, cfg.rna_dna_norm_witten_bell_output_files):
            self.assertTrue(utils.files_almost_equal(output_basename+file_extension, expected_output_file))

        # remove the temp file
        utils.remove_temp_folder(tempdir)
        
    def test_humann2_rna_dna_norm_log_tsv(self):
        """
        Test norm the tsv file entries from dna and rna input files with humann2_rna_dna_norm_table
        Test with log transform
        """
        
        # create a temp folder
        tempdir=utils.create_temp_folder("rna_dna_norm_log")
        output_basename=os.path.join(tempdir,"rna_dna_norm")
        
        # run the command
        utils.run_command(["humann2_rna_dna_norm","--input_dna",cfg.rna_dna_norm_dna_input,
                           "--input_rna",cfg.rna_dna_norm_rna_input,"--output_basename",
                           output_basename,"--log_transform"])
        
        # check the output files are as expected
        for file_extension, expected_output_file in zip(cfg.rna_dna_norm_file_names, cfg.rna_dna_norm_log_output_files):
            self.assertTrue(utils.files_almost_equal(output_basename+file_extension, expected_output_file))

        # remove the temp file
        utils.remove_temp_folder(tempdir)
        
    def test_humann2_rna_dna_norm_log_10_tsv(self):
        """
        Test norm the tsv file entries from dna and rna input files with humann2_rna_dna_norm_table
        Test with log transform with base 10
        """
        
        # create a temp folder
        tempdir=utils.create_temp_folder("rna_dna_norm_log_10")
        output_basename=os.path.join(tempdir,"rna_dna_norm")
        
        # run the command
        utils.run_command(["humann2_rna_dna_norm","--input_dna",cfg.rna_dna_norm_dna_input,
                           "--input_rna",cfg.rna_dna_norm_rna_input,"--output_basename",
                           output_basename,"--log_transform", "--log_base","10"])
        
        # check the output files are as expected
        # allow for varying precision in the calculations with almost equal
        for file_extension, expected_output_file in zip(cfg.rna_dna_norm_file_names, cfg.rna_dna_norm_log_10_output_files):
            self.assertTrue(utils.files_almost_equal(output_basename+file_extension, expected_output_file))

        # remove the temp file
        utils.remove_temp_folder(tempdir)
        
    def test_humann2_strain_profile_tsv(self):
        """
        Test the tsv file entries running humann2_strain_profile
        Test with critical mean and critical count values
        """
        
        # create a temp folder
        tempdir=utils.create_temp_folder("strain_profile")
        
        # move to this folder as the output files will be created in the current working folder
        current_working_directory=os.getcwd()
        try:
            os.chdir(tempdir)
        except EnvironmentError:
            print("Warning: Unable to move to temp directory: " + tempdir)
        
        # run the command
        utils.run_command(["humann2_strain_profiler","--input",cfg.strain_profile_input,
                           "--critical_mean","1","--critical_count","2"])
        
        # check the output files are as expected
        # allow for varying precision in the calculations with almost equal
        for file, expected_output_file in zip(cfg.strain_profile_file_names, cfg.strain_profile_m1_n2_output_files):
            self.assertTrue(utils.files_almost_equal(os.path.join(tempdir,file), expected_output_file))

        # return to original working directory
        os.chdir(current_working_directory)

        # remove the temp file
        utils.remove_temp_folder(tempdir)
        
    def test_humann2_merge_abundance_tsv(self):
        """
        Test the tsv gene families and pathway abundance file entries with humann2_merge_abundance_tables
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_merge_abundance_tables","--input-genes",cfg.merge_abundance_genefamilies_input,
                           "--input-pathways",cfg.merge_abundance_pathways_input,"--output",
                           new_file])
        
        # check the output file is as expected
        # allow for varying precision in the calculations with almost equal
        self.assertTrue(utils.files_almost_equal(new_file, cfg.merge_abundance_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_merge_abundance_remove_taxonomy_tsv(self):
        """
        Test the tsv gene families and pathway abundance file entries with humann2_merge_abundance_tables
        Test with the remove taxonomy option which stratifies by pathway then gene instead of
        stratifying by pathway, taxonomy, then gene
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_merge_abundance_tables","--input-genes",cfg.merge_abundance_genefamilies_input,
                           "--input-pathways",cfg.merge_abundance_pathways_input,"--output",
                           new_file,"--remove-taxonomy"])
        
        # check the output file is as expected
        # allow for varying precision in the calculations with almost equal
        self.assertTrue(utils.files_almost_equal(new_file, cfg.merge_abundance_remove_taxonomy_output))

        # remove the temp file
        utils.remove_temp_file(new_file)
