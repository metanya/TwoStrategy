import os
from collections import Counter

from Bio import Entrez, SeqIO

from record.helpers import get_interregions
from .constants import ENTREZ_EMAIL, GENE_BANK_FOLDER, GENES


def assert_gb_folder():
    if not os.path.exists(GENE_BANK_FOLDER):
        os.makedirs(GENE_BANK_FOLDER)


def assert_sum_of_genes(dictionary):
    gene_num = dictionary["gene"]
    number_of_occurences = sum(dictionary[gene_type] for gene_type in GENES)
    assert gene_num == number_of_occurences, "gene number is: {0} and number of the genes " \
                                             "types is: {1}".format(gene_num, number_of_occurences)


def assert_percentage(number):
    assert 0 <= number <= 100


def get_protein_gc_number(seq):
    counter_g = seq.count('G')
    counter_c = seq.count('C')
    return counter_c + counter_g


def assert_number_of_intergenes_are_greater_than_genes(size_of_genes, size_of_intergene):
    assert size_of_genes <= size_of_intergene


class Record:

    def __init__(self, record_id,record_family, parser):
    #def __init__(self, record_id, record_family, parser):
        self.main_attributes_dictionary = {}
        self.record_id = record_id
        self.record_family = record_family
        self.create_genbank_file()
        record_content = self.get_record_content()
        self.taxonomy = record_content.annotations['taxonomy']
        self.df = parser.get_data_frame('data\\csv\\{}.csv'.format(record_id), record_content)
        self.get_main_attributes(record_content)

    def search(self):
        Entrez.email = ENTREZ_EMAIL
        handle = Entrez.esearch(db="nucleotide", term=self.record_id)
        record = Entrez.read(handle)
        a = record["IdList"][0]

        record = Entrez.read(Entrez.elink(dbfrom="nucleotide", id=a))
        print(record)

    def get_genbank_record(self):
        with Entrez.efetch(db="nucleotide", id=self.record_id, rettype="gb", retmode="full",
                           usehistory="true", style='gbwithparts') as handle:
            list_of_records = []
            for record in SeqIO.parse(handle, "genbank"):
                list_of_records.append(record)
                print()
            return list_of_records[0]

    def get_record_content(self):
        assert os.path.exists(GENE_BANK_FOLDER + '{}.gb'.format(self.record_id))
        file_name = GENE_BANK_FOLDER + '{}.gb'.format(self.record_id)
        with open(file_name, "r") as handle:
            for i, record_gb in enumerate(SeqIO.parse(handle, "genbank")):
                return record_gb  # next(record_gb) # the last record

    def create_genbank_file(self):
        assert_gb_folder()
        # pubDateEnd = "2012/12/27"
        # pubDateStart = "2003/7/25"
        # searchTerm = f'("{pubDateStart}"[Publication Date]: "{pubDateEnd}"[Publication Date])'
        if not os.path.exists(GENE_BANK_FOLDER + '{}.gb'.format(self.record_id)):  # if the file not exists
            Entrez.email = ENTREZ_EMAIL

            with Entrez.efetch(db="nucleotide", id=self.record_id, rettype="gbwithparts",
                               retmode="text") as handle:  # ,  term=searchTerm
                with open(GENE_BANK_FOLDER + '{}.gb'.format(self.record_id), "w") as out_handle:
                    out_handle.write(handle.read())
                print("The file: {}.gb created".format(self.record_id))

    def get_main_attributes(self, record_content):
        self.main_attributes_dictionary["name"]=self.record_id ##
        self.main_attributes_dictionary["family"] = self.record_family###

        genome_size = self.df["length"][0]
        self.main_attributes_dictionary["genome_size"] = genome_size
        genes_counter_dictionary = Counter(self.df["type"])


        # assert_sum_of_genes(genes_counter_dictionary)
        self.main_attributes_dictionary.update(genes_counter_dictionary)
        genes_seq_len = self.df.loc[self.df['type'] == 'gene', 'length'].sum()
        self.main_attributes_dictionary["gene_length_in_genome"] = genes_seq_len
        gene_length_percent_in_genome = (genes_seq_len / (genome_size * 2)) * 100
        assert_percentage(gene_length_percent_in_genome)
        self.main_attributes_dictionary["%_gene_length_in_genome"] = gene_length_percent_in_genome

        assert_percentage(self.df["gc_percentage"][0])
        self.main_attributes_dictionary["%_GC_in_genome"] = self.df["gc_percentage"][0]  # to find

        GC_in_genes_number = self.df.loc[self.df['type'] == 'gene', 'gc_number'].sum()
        percentage_of_GC_in_genes = (GC_in_genes_number / genes_seq_len) * 100
        assert_percentage(percentage_of_GC_in_genes)
        self.main_attributes_dictionary["%_GC_in_genes"] = percentage_of_GC_in_genes

        intergenes, length_of_intergenes = get_interregions(record_content)
        # assert_number_of_intergenes_are_greater_than_genes(genes_counter_dictionary['gene'], len(intergenes))

        self.main_attributes_dictionary["intergene_length_in_genome"] = length_of_intergenes
        self.main_attributes_dictionary["%_intergene_length_in_genome"] = (length_of_intergenes / (
                genome_size * 2)) * 100
        length_of_intergenes_gc = sum(get_protein_gc_number(intergene.seq.upper()) for intergene in intergenes)
        self.main_attributes_dictionary["%_GC_in_intergene"] = (length_of_intergenes_gc
                                                                   / length_of_intergenes) * 100
