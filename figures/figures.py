import os
import numpy
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn
from scipy.stats.stats import pearsonr
from collections import OrderedDict
import Bio.Data.CodonTable as Codon
import pandas as pds
from figures.constants import MYOVIRIDAE, PODOVIRIDAE, SYNECHOCOCCUS, PROCHLOROCOCCUS, HL_PROCHLOROCOCCUS, LL_PROCHLOROCOCCUS
from figures.constants import BLUE, RED, BLACK, NONE
from figures.constants import T_RNA, PERCENTAGE_GC_IN_GENES, FAMILY, NAME

def show_relative_position(data, attributes, record_family, record_id):
    trna_data = data[data['type'] == 'tRNA']
    if record_family == "Myoviridae":
        seaborn.stripplot(x="relative_position", data=trna_data, dodge=True, palette=["red"], marker="*",
                          linewidth=1)
    elif record_family == "Podoviridae":
        seaborn.stripplot(x="relative_position", data=trna_data, dodge=True, palette=["blue"], marker="*",
                          linewidth=1)
    elif record_family == "LL-Prochlorococcus" or record_family == "HL-Prochlorococcus":
        seaborn.stripplot(x="relative_position", data=trna_data, dodge=True, palette=["white"], marker="^",
                          linewidth=1)
    else:  # Synechococcus
        seaborn.stripplot(x="relative_position", data=trna_data, dodge=True, palette=["white"], marker="o",
                          linewidth=1)
    plt.xlim([0, 1])
    # attributes.loc[attributes['name'] == record_id, 'percentage_GC_in_genome'][0]
    # a=attributes.loc[attributes['name'] == str(record_id), 'percentage_GC_in_genes'].tolist()
    plt.title(record_id + " | (%GC in genome:" + str(
        round(attributes.loc[attributes['name'] == record_id, 'percentage_GC_in_genome'][0],
              3)) + ") | (" + record_family + ")")
    # attributes.loc[attributes['name'] == record_id, 'percentage_GC_in_genes']

    # plt.title(record_id +" ("+record_family+")")

    if not os.path.exists("data\\png\\relative_position\\"):
        os.makedirs("data\\png\\relative_position\\")
    plt.savefig('data\\png\\relative_position\\{}.png'.format(record_id))

    plt.show()


def get_vals_by_family_and_attribute(main_attributes: pd.DataFrame, family: str, column: str) -> list[int]:
    return main_attributes.loc[main_attributes[FAMILY] == family, column].tolist()


def add_attributes_to_lists(host_list, virus_infect_host_gc_percent, virus_infect_host_trna, row, host_family: str, attribute1: str,
      attribute2):
    if any(host_family in s for s in host_list):
        virus_infect_host_gc_percent.append(row[attribute1])
        if not isinstance(row[attribute2], int):
            row[attribute2] = 0
        virus_infect_host_trna.append(row[attribute2])


class Figures:
    def __init__(self):
        pass
        # self.records = records
        # self.types = types
        # self.viruses_and_hosts_they_infect = viruses_and_hosts_they_infect
        # self.main_attributes_all_species = main_attributes_all_species

    @staticmethod
    def bar_chart_histogram(data: dict, ylabel, title, file_name, colors=['blue', 'red', 'white', 'white']):
        keys = data.keys()
        x_pos = np.arange(len(keys))  # number of bars
        means = [x[0] for x in data.values()]  # means or heights of the bars
        errors = [x[1] for x in data.values()]  # the error bars and the standard deviations

        # Build the plot
        fig, ax = plt.subplots()
        ax.bar(x_pos, means, yerr=errors, align='center', alpha=0.5, ecolor='black', capsize=10, color=colors,
               edgecolor="black")
        ax.set_ylabel(ylabel)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(keys)
        ax.set_title(title)
        ax.yaxis.grid(True)

        # Save the figure and show
        plt.tight_layout()
        plt.savefig(file_name)
        plt.show()

    def get_mean_and_std(self, records, types):
        mean_and_std_of_types = {}
        for organism_type in types:
            mean_and_std_of_types[organism_type] = self.get_mean_and_std_by_type(records, organism_type)
        return mean_and_std_of_types

    @staticmethod
    def get_mean_and_std_by_type(records, organism_type):
        sizes = numpy.array([])
        for record in records:
            if organism_type in record.taxonomy:
                sizes = numpy.append(sizes, record.main_attributes['genome_size'])
        # assert len(sizes) == 0, "sizes list is empty"
        if len(sizes) != 0:
            return sizes.mean(), sizes.std()

    @staticmethod
    def scatter_plot(main_attributes, viruses_and_hosts_they_infect):
        # Synechococcus organisms:
        syn_gc = get_vals_by_family_and_attribute(main_attributes, SYNECHOCOCCUS, PERCENTAGE_GC_IN_GENES)
        syn_trna = get_vals_by_family_and_attribute(main_attributes, SYNECHOCOCCUS, T_RNA)
        # LL-Prochlorococcus organisms:
        llproch_gc = get_vals_by_family_and_attribute(main_attributes, LL_PROCHLOROCOCCUS,PERCENTAGE_GC_IN_GENES)
        llproch_trna = get_vals_by_family_and_attribute(main_attributes, LL_PROCHLOROCOCCUS, T_RNA)
        # HL-Prochlorococcus organisms:
        hlproch_gc = get_vals_by_family_and_attribute(main_attributes, HL_PROCHLOROCOCCUS,PERCENTAGE_GC_IN_GENES)
        hlproch_trna = get_vals_by_family_and_attribute(main_attributes, HL_PROCHLOROCOCCUS, T_RNA)

        # Myoviridae organisms:
        myo_df = main_attributes.loc[main_attributes[FAMILY] == MYOVIRIDAE]
        # Myoviridae (host: Synechococcus):
        myo_infect_syn_gc, myo_infect_syn_trna = [], []
        # Myoviridae (host: Prochlorococcus):
        myo_infect_llproch_gc, myo_infect_llproch_trna = [], []
        myo_infect_hlproch_gc, myo_infect_hlproch_trna = [], []
        myo_infect_proch_gc, myo_infect_proch_trna = [], []

        for index, row in myo_df.iterrows():
            myo_host_list = viruses_and_hosts_they_infect.get(row[NAME])
            add_attributes_to_lists(myo_host_list, myo_infect_syn_gc, myo_infect_syn_trna, row, SYNECHOCOCCUS, PERCENTAGE_GC_IN_GENES, T_RNA)
            add_attributes_to_lists(myo_host_list, myo_infect_llproch_gc, myo_infect_llproch_trna, row, LL_PROCHLOROCOCCUS, PERCENTAGE_GC_IN_GENES, T_RNA)
            add_attributes_to_lists(myo_host_list, myo_infect_hlproch_gc, myo_infect_hlproch_trna, row, HL_PROCHLOROCOCCUS, PERCENTAGE_GC_IN_GENES, T_RNA)
            add_attributes_to_lists(myo_host_list, myo_infect_proch_gc, myo_infect_proch_trna, row, PROCHLOROCOCCUS, PERCENTAGE_GC_IN_GENES, T_RNA)

        # Podoviridae organisms:
        podo_df = main_attributes.loc[main_attributes[FAMILY] == PODOVIRIDAE]
        # Podoviridae (host: Synechococcus):
        podo_infect_syn_gc, podo_infect_syn_trna = [], []
        # Podoviridae (host: Prochlorococcus):
        podo_infect_llproch_gc, podo_infect_llproch_trna = [], []
        podo_infect_hlproch_gc, podo_infect_hlproch_trna = [], []

        for index, row in podo_df.iterrows():
            podo_host_list = viruses_and_hosts_they_infect.get(row[NAME])
            add_attributes_to_lists(podo_host_list, podo_infect_syn_gc, podo_infect_syn_trna, row, SYNECHOCOCCUS, PERCENTAGE_GC_IN_GENES, T_RNA)
            add_attributes_to_lists(podo_host_list, podo_infect_llproch_gc, podo_infect_llproch_trna, row, LL_PROCHLOROCOCCUS, PERCENTAGE_GC_IN_GENES, T_RNA)
            add_attributes_to_lists(podo_host_list, podo_infect_hlproch_gc, podo_infect_hlproch_trna, row, HL_PROCHLOROCOCCUS, PERCENTAGE_GC_IN_GENES, T_RNA)

        plt.scatter(syn_gc, syn_trna, c=NONE, marker="o", edgecolor=BLACK)

        plt.scatter(llproch_gc, llproch_trna, c=NONE, marker="^", edgecolor=BLACK)
        plt.scatter(hlproch_gc, hlproch_trna, c=NONE, marker="^", edgecolor=BLACK)

        plt.scatter(myo_infect_syn_gc, myo_infect_syn_trna, c=RED, linewidths=2, marker="o", edgecolor=NONE)
        plt.scatter(myo_infect_llproch_gc, myo_infect_llproch_trna, c=RED, linewidths=2, marker="^", edgecolor=NONE)
        plt.scatter(myo_infect_hlproch_gc, myo_infect_hlproch_trna, c=RED, linewidths=2, marker="^", edgecolor=NONE)
        plt.scatter(myo_infect_proch_gc, myo_infect_proch_trna, c=RED, linewidths=2, marker="^", edgecolor=NONE)

        plt.scatter(podo_infect_syn_gc, podo_infect_syn_trna, c=BLUE, marker="o", edgecolor=NONE)
        plt.scatter(podo_infect_llproch_gc, podo_infect_llproch_trna, c=BLUE, marker="^", edgecolor=NONE)
        plt.scatter(podo_infect_hlproch_gc, podo_infect_hlproch_trna, c=BLUE, marker="^", edgecolor=NONE)

        plt.title("(A) The tRNA Gene Copy Number of various cyanobacteria and cyanophages against their GC Content")
        plt.legend(
            [SYNECHOCOCCUS, "LL-Prochlorococcus", "HL-Prochlorococcus", ("Myo (host: %s)" % SYNECHOCOCCUS),
             "Myo (host: LL-Prochlorococcus)", "Myo (host: HL-Prochlorococcus)", "Myo (host: Prochlorococcus)",
             ("Podo (host: %s)" % SYNECHOCOCCUS), "Podo (host: LL-Prochlorococcus)", "Podo (host: HL-Prochlorococcus)"])
        plt.xlabel("GC content (%)")
        plt.ylabel("tRNA Gene Copy Number")
        plt.show()

    @staticmethod
    def get_codon_table(table_id: int):
        return Codon.generic_by_id[table_id]

    @staticmethod
    def get_start_codon(seq_part, codon_table):
        start_codon = -1
        if seq_part[0:3] in codon_table.start_codons:
            start_codon = 0
        elif seq_part[1:4] in codon_table.start_codons:
            start_codon = 1
        elif seq_part[2:5] in codon_table.start_codons:
            start_codon = 2
        print(f'Start codon is: {start_codon}')
        return start_codon

    @staticmethod
    def get_frequency_of_codons(cds_list, start_codon_list):
        print('Get codons frequency')
        frequency = OrderedDict()
        for index, cds in enumerate(cds_list):
            if start_codon_list[index] is not None:  # when abnormal
                for iteration in range(start_codon_list[index], len(cds), 3):
                    codon = cds[iteration: iteration + 3]
                    if codon in frequency:
                        frequency[codon] += 1
                    else:
                        if len(codon) == 3:
                            frequency[codon] = 1
        return OrderedDict(sorted(frequency.items()))

    # @staticmethod
    # def get_correlation(list_one, list_two):
    #     return pearsonr(list_one,list_two)[0]

    @staticmethod
    def get_correlation(dict_one: {}, dict_two: {}):
        return pearsonr([f[1] for f in list(dict_one.items())], [f[1] for f in list(dict_two.items())])[0]

        # # seaborn.stripplot( x="Correlation", y="between", hue=None, data=None, order=None, hue_order=None,
        # jitter=True, dodge=True,
        # #                   orient=None, color=None, palette=None, size=5, edgecolor='gray', linewidth=0, ax=None)

        # data = {"species": ["Myo-Myo", "Myo-Myo", "Myo-Myo", "Myo-Podo", "HLProch-Syne"],
        #                     "Correlation": [0.5, 0.4, 0, 0.3, 0.6]};
        # df = pds.DataFrame(data=data);
        # print(df);
        # seaborn.stripplot(x="Correlation", y="species", data=df);
        # plt.show();

    # def stripchart(self, frequencies_dict, main_attributes, stop_codons):
    #     species = []
    #     frequencies = []
    #     species_vs_species = []
    #     correlation = []
    #
    #     for stop_codon in stop_codons:
    #         for frequency_dict in frequencies_dict.values():
    #             frequency_dict.pop(stop_codon, None)
    #
    #     for key, value in frequencies_dict.items():
    #         species.append(main_attributes.loc[main_attributes['name'] == key, 'family'][0])
    #         frequencies.append(value)
    #
    #     for i in range(len(species) - 1):
    #         for j in range(i + 1, len(species)):
    #             species_vs_species.append(str(species[i]) + "+" + str(species[j]))
    #             correlation.append(self.get_correlation(frequencies[i], frequencies[j]))
    #
    #     data = {"species": species_vs_species, "Correlation": correlation}
    #     df = pds.DataFrame(data=data)
    #     seaborn.stripplot(x="Correlation", y="species", data=df)
    #     plt.grid(color='green', linestyle='--', linewidth=0.5)
    #     plt.show()
    # ---------------------------------------------------
    def draw_on_stripplot(self, data, color, marker):
        # plt.scatter(x="Correlation", y="species", data=data, c=color, marker=marker, edgecolor="black", linewidth=0.5,y_jitter=0.2)
        order = ["HL-Prochlorococcus+HL-Prochlorococcus", "LL-Prochlorococcus+LL-Prochlorococcus",
                 "Synechococcus+Synechococcus", "Podoviridae+Podoviridae", "Myoviridae+Myoviridae",
                 "Pro-infecting Podoviridae+HL-Prochlorococcus", "Pro-infecting Podoviridae+LL-Prochlorococcus",
                 "Syn-infecting Podoviridae+Synechococcus", "Pro-infecting Myoviridae+HL-Prochlorococcus",
                 "Pro-infecting Myoviridae+LL-Prochlorococcus", "Syn-infecting Myoviridae+Synechococcus"]
        seaborn.stripplot(x="Correlation", y="species", data=data, order=order, dodge=True, palette=color,
                          marker=marker,
                          linewidth=0.5, jitter=0.35, edgecolor='black')

    def add_data_to_stripplot(self, df, same_species):
        if same_species:
            self.draw_on_stripplot(df[(df['species'] == "HL-Prochlorococcus+HL-Prochlorococcus")], ["white"], "^")
            self.draw_on_stripplot(df[(df['species'] == "LL-Prochlorococcus+LL-Prochlorococcus")], ["white"], "^")
            self.draw_on_stripplot(df[(df['species'] == "Synechococcus+Synechococcus")], ["white"], "o")
            self.draw_on_stripplot(df[(df['species'] == "Podoviridae+Podoviridae")], ["blue"], "*")
            self.draw_on_stripplot(df[(df['species'] == "Myoviridae+Myoviridae")], ["red"], "*")
        else:
            self.draw_on_stripplot(df[(df['species'] == "Pro-infecting Podoviridae+HL-Prochlorococcus")], ["blue"], "^")
            self.draw_on_stripplot(df[(df['species'] == "Pro-infecting Podoviridae+LL-Prochlorococcus")], ["blue"], "^")
            self.draw_on_stripplot(df[(df['species'] == "Syn-infecting Podoviridae+Synechococcus")], ["blue"], "o")
            self.draw_on_stripplot(df[(df['species'] == "Pro-infecting Myoviridae+HL-Prochlorococcus")], ["red"], "^")
            self.draw_on_stripplot(df[(df['species'] == "Pro-infecting Myoviridae+LL-Prochlorococcus")], ["red"], "^")
            self.draw_on_stripplot(df[(df['species'] == "Syn-infecting Myoviridae+Synechococcus")], ["red"], "o")

    def test_stripchart(self, pairs_list, frequencies_dict, main_attributes,
                        viruses_and_hosts_they_infect, stop_codons):
        species = []
        frequencies = []
        infecting = []

        # remove the stop codons
        for stop_codon in stop_codons:
            for frequency_dict in frequencies_dict.values():
                frequency_dict.pop(stop_codon, None)

        # remove codons that contain 'N'
        for frequency_dict in frequencies_dict.values():
            for key in frequency_dict.copy().keys():
                if 'N' in key:
                    frequency_dict.pop(key, None)

        for key, value in frequencies_dict.items():
            species.append(main_attributes.loc[main_attributes['name'] == key, 'family'][0])
            frequencies.append(value)
            if species[-1] in ["Myoviridae", "Podoviridae"]:  # if the item is virus
                infecting.append(viruses_and_hosts_they_infect[key])  # add the hosts it infect
            else:
                infecting.append([""])

        df_same_species_pairs = self.get_df_same_species_pairs(species, frequencies, pairs_list)

        df_virus_host_pairs = self.get_df_virus_host_pairs(species, frequencies, pairs_list, infecting)

        self.add_data_to_stripplot(df_same_species_pairs, True)
        self.add_data_to_stripplot(df_virus_host_pairs, False)

        plt.grid(color='green', linestyle='--', linewidth=0.5)
        plt.xlabel("Correlation")
        plt.ylabel("species Vs. species")
        plt.show()

    def get_df_same_species_pairs(self, species, frequencies, pairs_list):
        species_vs_species = []
        correlation = []
        for i in range(len(species) - 1):
            for j in range(i + 1, len(species)):
                species_vs_species.append(str(species[i]) + "+" + str(species[j]))
                correlation.append(self.get_correlation(frequencies[i], frequencies[j]))

        data = {"species": species_vs_species, "Correlation": correlation}
        df = pds.DataFrame(data=data)
        return df[df['species'].isin(pairs_list)]

    def get_df_virus_host_pairs(self, species, frequencies, pairs_list, infecting):
        species_vs_species = []
        correlation = []

        for i in range(len(species) - 1):
            for j in range(i + 1, len(species)):
                if species[i] in ["Myoviridae", "Podoviridae"]:  # if the item is virus
                    if ("LL-Prochlorococcus" in infecting[i]) | ("HL-Prochlorococcus" in infecting[i]) | (
                            "Prochlorococcus" in infecting[i]):
                        species_vs_species.append("Pro-infecting " + str(species[i]) + "+" + str(species[j]))
                        correlation.append(self.get_correlation(frequencies[i], frequencies[j]))
                    if "Synechococcus" in infecting[i]:
                        species_vs_species.append("Syn-infecting " + str(species[i]) + "+" + str(species[j]))
                        correlation.append(self.get_correlation(frequencies[i], frequencies[j]))
                else:
                    species_vs_species.append(str(species[i]) + "+" + str(species[j]))
                    correlation.append(self.get_correlation(frequencies[i], frequencies[j]))

        data = {"species": species_vs_species, "Correlation": correlation}
        return pds.DataFrame(data=data)

    def test_show_stripchart(self):
        data1 = {"days": ["Monday", "Monday"], "Pay": [6, 9]}
        data2 = {"days": ["Thursday", "Thursday"], "Pay": [5, 4]}

        all_days = ["Monday", "Thursday"]

        seaborn.stripplot(x="Pay", y="days", data=data1, order=all_days, dodge=True, palette=["red"], marker="^",
                          linewidth=1)
        seaborn.stripplot(x="Pay", y="days", data=data2, order=all_days, dodge=True, palette=["blue"], marker="o",
                          linewidth=1)

        plt.grid(color='green', linestyle='--', linewidth=0.5)
        plt.tight_layout()
        plt.show()

        # data1 = {"days": ["Monday", "Monday"], "Pay": [6, 9]}
        # data2 = {"days": ["Thursday", "Thursday"], "Pay": [5, 4]}
        #
        # plt.scatter(x="Pay", y="days", data=data1, c="red", marker="^", linewidth=1)
        # plt.scatter(x="Pay", y="days", data=data2, c="blue", marker="o", linewidth=1)
        #
        # plt.grid(color='green', linestyle='--', linewidth=0.5)
        # plt.ylim(len(plt.yticks()) - 0.5, -0.5)
        # plt.tight_layout()
        # plt.show()
