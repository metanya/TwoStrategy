from Bio import Entrez, SeqIO
from .constants import ENTREZ_EMAIL, GENE_BANK_FOLDER
import os


def assert_gb_folder():
    if not os.path.exists(GENE_BANK_FOLDER):
        os.makedirs(GENE_BANK_FOLDER)


class Record:
    def __init__(self, record_id):
        self.record_id = record_id

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
        file_name = GENE_BANK_FOLDER + '{}.gb'.format(self.record_id)
        with open(file_name, "r") as handle:
            for i, record_gb in enumerate(SeqIO.parse(handle, "genbank")):
                print('Record number: {}\n============='.format(i))
                print(record_gb)
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
