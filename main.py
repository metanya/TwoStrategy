from figures.figures import Figures
from records.records import Records


def species_names_dictionary_old():
    return {
        # "KM034562.1":"KM034562.1",
        # "NC_045512.2": "NC_045512.2", # cororna virus
        # '''***Cyanophages:***'''
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

        # '''***Cyanobacteria:***'''
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
        # "MIT_9211": "NC_009976.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        "MIT_9303": "NC_008820.1",  # circular CON 25-NOV-2016 !!!!!!!!!!!!!
        # "MIT_9313": "NC_005071.1",  # circular CON 13-DEC-2020 !!!!!!!!!!!!!

        # ------Synechococcus:------
        "WH_8102": "NC_005070.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        "WH_7803": "NC_009481.1",  # circular CON 17-APR-2017 !!!!!!!!!!!!!

    }


def species_names_dictionary():
    # species_names_dictionary = {{'virus1':'virus1 family'},{'virus2':'virus1 family'},....}
    return {
        # "KM034562.1":"KM034562.1",
        # "NC_045512.2": "NC_045512.2", # cororna virus
        # '''***Cyanophages:***'''
        # ------Moyviruses:------
        "NC_006884.1": "Myoviridae",  # [Phage: Myoviridae P-SSM4], circular PHG 26-MAR-2010
        "NC_006883.1": "Myoviridae",  # [Phage: Myoviridae P_SSM2], circular PHG 26-MAR-2010
        "NC_008296.1": "Myoviridae",  # [Phage: Myoviridae Syn9], circular PHG 28-NOV-2007
        "FM207411.1": "Myoviridae",  # [Phage: Myoviridae S_RSM4], circular PHG 22-SEP-2009
        "NC_006820.1": "Myoviridae",  # [Phage: Myoviridae S_PM2], circular PHG 11-OCT-2021 !
        # # ------Podoviruses:------
        "NC_006882.1": "Podoviridae",  # [Phage: Podoviridae P_SSP7], circular PHG 19-NOV-2010
        "NC_009531.1": "Podoviridae",  # [Phage: Podoviridae Syn5], linear   PHG 20-DEC-2020 !
        "NC_003390.1": "Podoviridae",  # [Phage: Podoviridae P60], linear   PHG 17-APR-2009

        # '''***Cyanobacteria:***'''
        # ------Prochlorococcus:------
        # - - - - - - HL-Prochlorococcus:
        "NC_005072.1": "HL-Prochlorococcus",  # [Cyanobacteria: HL-Prochlorococcus MED4], circular CON 10-OCT-2021  !
        "NC_008817.1": "HL-Prochlorococcus",  # [Cyanobacteria: HL-Prochlorococcus MIT_9515], circular CON 10-OCT-2021 !
        "NC_007577.1": "HL-Prochlorococcus",  # [Cyanobacteria: HL-Prochlorococcus MIT_9312], circular CON 10-OCT-2021 !
        "NC_009840.1": "HL-Prochlorococcus",  # [Cyanobacteria: HL-Prochlorococcus MIT_9215], circular CON 10-OCT-2021 !
        # - - - - - - LL-Prochlorococcus:
        "NC_008819.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus NATL1A], circular CON 10-OCT-2021 !
        "NC_007335.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus NATL2A], circular BCT 29-NOV-2007
        "NC_005042.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus SS120], circular CON 10-OCT-2021 !
        # "NC_009976.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus MIT_9211], circular CON 10-OCT-2021 !
        "NC_008820.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus MIT_9303], circular CON 25-NOV-2016 !
        # "NC_005071.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus MIT_9313], circular CON 13-DEC-2020 !

        # ------Synechococcus:------
        "NC_005070.1": "Synechococcus",  # [Cyanobacteria: Synechococcus WH_8102], circular CON 10-OCT-2021 !
        "NC_009481.1": "Synechococcus",  # [Cyanobacteria: Synechococcus WH_7803], circular CON 17-APR-2017 !

    }


def viruses_and_hosts_they_infect():
    # viruses_and_hosts_they_infect = {'virus1':['infected1','infected2'],'virus2':['infected1','infected2'],....}
    return {
        # '''***Cyanophages:***'''
        # ------Moyviruses:------
        "NC_006884.1": ['LL-Prochlorococcus', 'HL-Prochlorococcus'],
        # [Phage: Myoviridae P-SSM4], circular PHG 26-MAR-2010
        "NC_006883.1": ['LL-Prochlorococcus'],  # [Phage: Myoviridae P_SSM2], circular PHG 26-MAR-2010
        "NC_008296.1": ['Synechococcus', 'LL-Prochlorococcus', 'HL-Prochlorococcus'],
        # [Phage: Myoviridae Syn9], circular PHG 28-NOV-2007
        "FM207411.1": ['Synechococcus'],  # [Phage: Myoviridae S_RSM4], circular PHG 22-SEP-2009
        "NC_006820.1": ['Synechococcus2'],  # [Phage: Myoviridae S_PM2], circular PHG 11-OCT-2021 !
        # ------Podoviruses:------
        "NC_006882.1": ['HL-Prochlorococcus'],  # [Phage: Podoviridae P_SSP7], circular PHG 19-NOV-2010
        "NC_009531.1": ['Synechococcus'],  # [Phage: Podoviridae Syn5], linear   PHG 20-DEC-2020 !
        "NC_003390.1": ['Synechococcus'],  # [Phage: Podoviridae P60], linear   PHG 17-APR-2009
    }


def records(dictionary):
    records_obj = Records(dictionary)
    records = records_obj.records
    main_attributes_all_species = records_obj.df_main_attributes_for_all_species
    return records_obj, records, main_attributes_all_species


def manageFigures(records_obj, viruses_and_hosts_they_infect):
    frequencies = {}
    figures = Figures(records_obj, ["Podoviridae", "Myoviridae", "Prochlorococcus", "Synechococcus"],
                      viruses_and_hosts_they_infect)
    x = figures.get_codon_table(11)
    for record in records_obj.records:
        start_codon = figures.get_start_codon(record.record_content.seq, x)
        codons = figures.get_frequency_of_codons(start_codon, record.record_content.seq)
        # codons_vector = figures.get_vector_of_dictionary(codons)
        frequencies[record.record_id] = codons
        # !#frequencies[record.record_family] = codons
        # !#families.append(record.family)
    # values = list(frequencies.values())
    # #!#keys = list(frequencies.keys())
    # #!#for
    #     #!#for
    # cor = figures.get_correlation(values[0], values[1])
    # #!#specieVsspecie=keys[0]+"-"+keys[1]
    # print(cor)
    # mean_and_std_of_types = figures.get_mean_and_std()
    #
    figures.stripchart(frequencies)


# figures.bar_chart_histogram(mean_and_std_of_types)  # figure1, B
#  figures.scatter_plot()  # mean_and_std_of_types, viruses_and_hosts_they_infect) # figure1, A


if __name__ == '__main__':
    species_names_dictionary = species_names_dictionary()
    viruses_and_hosts_they_infect = viruses_and_hosts_they_infect()

    records_obj, records, main_attributes_all_species = records(species_names_dictionary)

    manageFigures(records_obj, viruses_and_hosts_they_infect)
