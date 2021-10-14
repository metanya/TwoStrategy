import pandas as pd


def get_main_attributes(record):
    '''
    type_col = pd.read_csv('data\\csv\\{}.csv'.format(record.record_id))
    gene_num = type_col["type"].to_list().count('gene')-1
    cds_num = type_col["type"].to_list().count('CDS')-1
    trna_num = type_col["type"].to_list().count('tRNA')
    '''
    record.main_attributes_dictionary["gene_num"] = record.df["type"].to_list().count('gene') - 1
    record.main_attributes_dictionary["cds_num"] = record.df["type"].to_list().count('CDS') - 1
    record.main_attributes_dictionary["trna_num"] = record.df["type"].to_list().count('tRNA')
    genome_size = record.df["length"][1]
    record.main_attributes_dictionary["genome_size"] = genome_size
    record.main_attributes_dictionary["percentage_of_genes_in_genome"] = (record.main_attributes_dictionary[
                                                                              "gene_num"] / genome_size) * 100
    record.main_attributes_dictionary["percentage_of_intergene_in_genome"] = ((genome_size -
                                                                               record.main_attributes_dictionary[
                                                                                   "gene_num"]) / genome_size) * 100

    record.main_attributes_dictionary["percentage_of_GC_in_genome"] = record.df["gc_percentage"][1]
    GC_in_genes_number =record.df.loc[record.df['type']=='gene','gc_number'].sum()
    genes_seq_len =record.df.loc[record.df['type']=='gene','length'].sum()
    GC_in_intergene_number =(record.df["gc_number"][1])-GC_in_genes_number
    intergene_seq_len =(record.df["length"][1])-genes_seq_len

    record.main_attributes_dictionary["percentage_of_GC_in_genes"] = (GC_in_genes_number / genes_seq_len) * 100
    record.main_attributes_dictionary["percentage_of_GC_in_intergene"] = (GC_in_intergene_number / intergene_seq_len) * 100
