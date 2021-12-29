import os

from record.helpers import create_data_dictionary
from parsers.helpers import transpose, convert_dictionary_to_csv_file


class Parser:
    def __init__(self):
        print()

    @staticmethod
    def get_record_data(csv_name, gene_bank_file: str):
        features = create_data_dictionary(gene_bank_file)
        df = transpose(features)
        if not os.path.exists('data\\csv\\'):
            os.makedirs('data\\csv\\')
        convert_dictionary_to_csv_file(df, csv_name)
        return df  # added

    @staticmethod
    def transpose(dictionary):
        return transpose(dictionary)

    @staticmethod
    def to_csv(data_frame, path):
        convert_dictionary_to_csv_file(data_frame, path)
