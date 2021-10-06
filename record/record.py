from Bio import Entrez, SeqIO
# from constants import ENTREZ_EMAIL
from .constants import ENTREZ_EMAIL
import os


class Record:
    def __init__(self, record_id):
        self.record_id = record_id

    def get_genbank_record(self):
        with Entrez.efetch(db="nucleotide", id=self.record_id, rettype="gb", retmode="text") as handle:
            # results = handle.read()
            # print(type(results))
            # return results
            list_of_records = []
            for record in SeqIO.parse(handle, "genbank"):
                list_of_records.append(record)
            # print(type(list_of_records[0]))
            return list_of_records[0]

    def get_record_content(self):
        file_name = 'data\\gb\\{}.gb'.format(self.record_id)
        with open(file_name, "r") as handle:
            for i, record_gb in enumerate(SeqIO.parse(handle, "genbank")):
                print('Record number: {}\n============='.format(i))
                print(record_gb)
            return record_gb  # next(record_gb) # the last record

    def create_genbank_file(self):
        if not os.path.exists('data\\gb\\{}.gb'.format(self.record_id)):  # if the file not exists

            Entrez.email = ENTREZ_EMAIL

            with Entrez.efetch(db="nucleotide", id=self.record_id, rettype="gb", retmode="text") as handle:
                with open('data\\gb\\{}.gb'.format(self.record_id), "w") as out_handle:
                    out_handle.write(handle.read())
                print("The file: {}.gb created".format(self.record_id))