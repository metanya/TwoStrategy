# from record.helpers import get_main_attributes
from record.record import Record
from parsers.parser import Parser
from records.records import Records

if __name__ == '__main__':
    # species_dictionary = {"species_name1":"Refseq1" , "species_name2":"Refseq2",.....}
    species_names_dictionary = {
        '''***Cyanophages:***'''
        # ------Moyviruses:------
        "P_SSM4": "NC_006884.1",  # P-SSM4: (Refseq=NC006884)
        "P_SSM2": "NC_006883",  # P-SSM2: (Refseq=NC006883)
        "Syn9": "NC_008296",  # Syn9: (Refseq=NC008296)
        "S_RSM4": "FM207411",  # S-RSM4: (Refseq=FM207411)
        "S_PM2": "NC_006820",  # S-PM2: (Refseq=NC006820)
        # ------Podoviruses:------
        "P_SSP7": "NC_006882",  # P-SSP7: (Refseq=NC006882)
        "Syn5": "NC_009531",  # Syn5: (Refseq=NC009531)
        "P60": "NC_003390",  # P60: (Refseq=NC003390)

        '''***Cyanobacteria:***'''
        # ------Prochlorococcus:------
        # - - - - - - HL-Prochlorococcus:
        "MED4": "NC_005072",  # Prochlorococcus marinus subsp pastoris str. CCMP1986 (MED4): (Refseq=NC_005072)
        "MIT_9515": "NC_008817",  # Prochlorococcus marinus str. MIT 9515: (Refseq=NC_008817)
        "MIT_9312": "NC_007577",  # Prochlorococcus marinus str. MIT 9312: (Refseq=NC_007577)
        "MIT_9215": "NC_009840",  # Prochlorococcus marinus str. MIT 9215: (Refseq=NC_009840)
        # - - - - - - LL-Prochlorococcus:
        "NATL1A": "NC_008819",  # Prochlorococcus marinus str. NATL1A: (Refseq=NC_008819)
        "NATL2A": "NC_007335",  # Prochlorococcus marinus str. NATL2A: (Refseq=NC_007335)
        "SS120": "NC_005042",  # Prochlorococcus marinus subsp marinus str. CCMP1375 (SS120): (Refseq=NC_005042)
        "MIT_9211": "NC_009976",  # Prochlorococcus marinus str. MIT 9211: (Refseq=NC_009976)
        "MIT_9303": "NC_008820",  # Prochlorococcus marinus str. MIT 9303: (Refseq=NC_008820)
        "MIT_9313": "NC_005071",  # Prochlorococcus marinus str. MIT 9313: (Refseq=NC_005071)

        # ------Synechococcus:------
        "WH_8102": "NC_005070",  # Synechococcus sp. WH 8102:(Refseq=NC_005070)
        "WH_7803": "NC_009481",  # Synechococcus sp. WH 7803:(Refseq=NC_009481)

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
