import os

from parsers.helpers import create_data_dictionary, transpose, convert_dictionary_to_csv_file


class Parser:
    def __init__(self):
        print()

    @staticmethod
    def get_data_frame(csv_name, gene_bank_file: str):
        features = create_data_dictionary(gene_bank_file)
        df = transpose(features)
        if not os.path.exists('data\\csv\\'):
            os.makedirs('data\\csv\\')
        convert_dictionary_to_csv_file(df, csv_name)
        return df#added
