from parsers.parser import Parser
from record.helpers import get_main_attributes
from record.record import Record
from records.helpers import create_main_attributes_dictionary_for_all_species
import pandas as pd


class Records:
    def __init__(self, records_ids):
        self.records = []

        parser = Parser()

        for record_id in records_ids:
            obj = Record(record_id)
            obj.create_genbank_file()
            self.records.append(obj)
            obj.df = parser.get_data_frame('data\\csv\\{}.csv'.format(obj.record_id), obj.get_record_content())  # added
            get_main_attributes(obj)
            # print(obj.main_attributes_dictionary)

        main_attributes_dictionary_for_all_species = create_main_attributes_dictionary_for_all_species(self.records)
        # create csv file of all the main attributes, of all the species:
        df_main_attributes_for_all_species = pd.DataFrame(main_attributes_dictionary_for_all_species)
        df_main_attributes_for_all_species.to_csv('data\\csv\\main_attributes_for_all_species.csv')
