import pandas as pd


def transpose(features):
    return pd.DataFrame(features, columns=features.keys()).transpose()


def convert_dictionary_to_csv_file(df, csv_name_path):
    # fileUtils.verify_file_is_closed(csv_name_path)
    df.to_csv(csv_name_path)
