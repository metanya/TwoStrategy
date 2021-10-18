# from record.helpers import get_main_attributes
from record.record import Record
from parsers.parser import Parser
from records.records import Records

if __name__ == '__main__':
    # species_dictionary = {"species_name1":"Refseq1" , "species_name2":"Refseq2",.....}
    species_names_dictionary = {
        '''***Cyanophages:***'''
        # ------Moyviruses:------
        "P_SSM4": "NC_006884.1",  # circular PHG 26-MAR-2010
        "P_SSM2": "NC_006883.1",  # circular PHG 26-MAR-2010
        "Syn9": "NC_008296.1",  # circular PHG 28-NOV-2007
        "S_RSM4": "FM207411.1",  # circular PHG 22-SEP-2009
        "S_PM2": "NC_006820.1",  # circular PHG 11-OCT-2021 !!!!!!!!!!!!!
        # ------Podoviruses:------
        "P_SSP7": "NC_006882.1",  # circular PHG 19-NOV-2010
        "Syn5": "NC_009531.1",  # linear   PHG 20-DEC-2020 !!!!!!!!!!!!!
        "P60": "NC_003390.1",  # linear   PHG 17-APR-2009

        '''***Cyanobacteria:***'''
        # ------Prochlorococcus:------
        # - - - - - - HL-Prochlorococcus:
        "MED4": "NC_005072.1",  # circular CON 10-OCT-2021  !!!!!!!!!!!!!
        "MIT_9515": "NC_008817.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        "MIT_9312": "NC_007577.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        "MIT_9215": "NC_009840.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        # - - - - - - LL-Prochlorococcus:
        "NATL1A": "NC_008819.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        "NATL2A": "NC_007335.1",  # circular BCT 29-NOV-2007
        "SS120": "NC_005042.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        "MIT_9211": "NC_009976.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        "MIT_9303": "NC_008820.1",  # circular CON 25-NOV-2016 !!!!!!!!!!!!!
        "MIT_9313": "NC_005071.1",  # circular CON 13-DEC-2020 !!!!!!!!!!!!!

        # ------Synechococcus:------
        "WH_8102": "NC_005070.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        "WH_7803": "NC_009481.1",  # circular CON 17-APR-2017 !!!!!!!!!!!!!

    }

    #parser = Parser()

    records = Records(species_names_dictionary.values()).records

    '''
    for species in species_names_dictionary:
        obj = Record(species_names_dictionary[species])
        obj.create_genbank_file()
        obj.df = parser.get_data_frame('data\\csv\\{}.csv'.format(obj.record_id), obj.get_record_content())  # added
        #get_main_attributes(obj)


        # parser.get_data_frame('data\\csv\\{}.csv'.format(obj.record_id), obj.get_record_content())
        # species_obj_dictionary["{0}".format(obj.record_id)] = obj.record_details
'''
    #print(get_main_attributes(Record("FM207411")))
