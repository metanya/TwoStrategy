from Bio import Entrez
from constants import ENTREZ_EMAIL


class Record:
    def __init__(self, record_id):
        self.record_id = record_id

    Entrez.email = ENTREZ_EMAIL

    def get_genbank_record(self):
        with Entrez.efetch(db="nucleotide", id=self.record_id, rettype="gb", retmode="text") as handle:
            results = handle.read()
            return results
