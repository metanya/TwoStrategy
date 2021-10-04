from record.record import Record
from parsers.parser import Parser

if __name__ == '__main__':
    record = Record("AP006627")
    record.create_genbank_file()
    AP006627_record_content = record.get_record_content()

    parser = Parser()
    parser.get_data_frame("test.csv", AP006627_record_content)
