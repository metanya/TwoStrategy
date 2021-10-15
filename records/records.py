from parsers.parser import Parser
# from record.helpers import get_main_attributes
from record.record import Record
from records.helpers import create_main_attributes_dictionary_for_all_species
import pandas as pd


class Records:
    def __init__(self, records_ids):
        self.records = []

        parser = Parser()

        for record_id in records_ids:
            print(record_id)
            record = Record(record_id, parser)
            self.records.append(record)

        main_attributes_dictionary_for_all_species = create_main_attributes_dictionary_for_all_species(self.records)
        df_main_attributes_for_all_species = pd.DataFrame(main_attributes_dictionary_for_all_species)
        df_main_attributes_for_all_species.to_csv('data\\csv\\main_attributes_for_all_species.csv')
