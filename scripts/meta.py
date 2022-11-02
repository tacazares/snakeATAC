"""meta.py

This file contains the class function for manipulating the experimental meta data

"""
import os
import pandas as pd


class MetaTable():
    """
    An object to hold meta data information
    """
    def __init__(self,
                 meta_path,
                 output_dir):
        self.df = pd.read_table(meta_path)

        self.output_dir = output_dir

        self.all_tech_reps = self.df["srr"].unique().tolist()

        self.bio_rep_list = self.df["srx"].unique().tolist()

        self.condition_list = self.df["condition"].unique().tolist()

        self.srx_2_condition_dict = pd.Series(self.df.condition.values,index=self.df.srx).to_dict()

        self.name_df = self.nameBuilderDF()
    
    #def getConditionReps(self, wildcards):
        #CONDITION_REP_LIST = self.df[self.df["condition"] == wildcards.condition]["srx"]

        #return list(map(lambda condition_rep_id: os.path.join(self.output_dir, condition, sample, "alignments", "star", sample ".bam"), CONDITION_REP_LIST))
    
    def nameBuilderDF(self):
        name_builder_df = self.df.copy()

        name_builder_df["read"] = [[1,2] for _ in range(len(name_builder_df))]

        name_df = name_builder_df.explode("read")

        return name_df

    def getReplicateFastq_pe1(self, wildcards):
        TECH_REP_LIST = self.df[self.df["srx"] == wildcards.sample]["srr"]

        return list(map(lambda tech_rep_id: os.path.join(self.output_dir, self.srx_2_condition_dict[wildcards.sample], "replicate_data", wildcards.sample, "fasterq_dump", tech_rep_id + "_1.fastq.gz"), TECH_REP_LIST))

    def getReplicateFastq_pe2(self, wildcards):
        TECH_REP_LIST = self.df[self.df["srx"] == wildcards.sample]["srr"]

        return list(map(lambda tech_rep_id: os.path.join(self.output_dir, self.srx_2_condition_dict[wildcards.sample], "replicate_data", wildcards.sample, "fasterq_dump", tech_rep_id + "_2.fastq.gz"), TECH_REP_LIST))
