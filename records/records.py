from parsers.parser import Parser
from record.record import Record
from records.helpers import create_main_attributes_dictionary_for_all_species


class Records:
    def __init__(self, records):
        self.records = []

        parser = Parser()

        for record_id in records:
            print(record_id)#record name
            print(records[record_id])#record family

            record = Record(record_id,records[record_id], parser)
            self.records.append(record)

        main_attributes_dictionary_for_all_species = create_main_attributes_dictionary_for_all_species(self.records)
        self.df_main_attributes_for_all_species = parser.transpose(main_attributes_dictionary_for_all_species)
        parser.to_csv(self.df_main_attributes_for_all_species, 'data\\csv\\main_attributes_for_all_species.csv')