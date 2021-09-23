from record.record import Record
from parsers.parser import Parser

if __name__ == '__main__':
    record = Record("AY851612").get_genbank_record()
    parser = Parser()
    parser.get_data_frame("test.csv", record)