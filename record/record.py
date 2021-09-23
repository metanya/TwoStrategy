from Bio import Entrez, SeqIO
from constants import ENTREZ_EMAIL


class Record:
    def __init__(self, record_id):
        self.record_id = record_id

    Entrez.email = ENTREZ_EMAIL

    def get_genbank_record(self):
        with Entrez.efetch(db="nucleotide", id=self.record_id, rettype="gb", retmode="text") as handle:
            list_of_records = []
            for record in SeqIO.parse(handle, "genbank"):
                list_of_records.append(record)
            print(type(list_of_records[0]))
            return list_of_records[0]
