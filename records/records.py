import logging
import os

import pandas as pd

from parsers.parser import Parser
from record.constants import TRNA_SCAN_CSV_FOLDER
from record.record import Record
from records.helpers import get_main_attributes

# def assert_compare_trna_num(trna_num_records, trna_num_scan):
#     assert trna_num_records == trna_num_scan


class Records:
    def __init__(self, records: dict[str, str]):
        self.records = []
        parser = Parser()

        for record_id in records:
            self.records.append(Record(record_id, records[record_id], parser))

        main_attributes = get_main_attributes(self.records)
        self.attributes = parser.transpose(main_attributes)
        parser.to_csv(self.attributes, 'data\\csv\\main_attributes_for_all_species.csv')

        # for record_id in records:
        #     if os.path.exists(TRNA_SCAN_CSV_FOLDER + "{}.csv".format(record_id)):
        #         trna_num_records=self.attributes[['tRNA']].loc[record_id]
        #         df = pd.read_csv(TRNA_SCAN_CSV_FOLDER + "{}.csv".format(record_id))
        #         trna_num_scan=df['tRNA Type'].sum()
        #         assert_compare_trna_num(trna_num_records, trna_num_scan)
