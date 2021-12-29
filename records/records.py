import logging
from parsers.parser import Parser
from record.record import Record
from records.helpers import get_main_attributes


class Records:
    def __init__(self, records):
        self.records = []
        parser = Parser()

        for record_id in records:
            self.records.append(Record(record_id, records[record_id], parser))

        main_attributes = get_main_attributes(self.records)
        self.attributes = parser.transpose(main_attributes)
        parser.to_csv(self.attributes, 'data\\csv\\main_attributes_for_all_species.csv')
