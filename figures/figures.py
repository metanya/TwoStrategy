import numpy
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import colors
from matplotlib.ticker import PercentFormatter


class Figures:
    def __init__(self,records_obj, types,viruses_and_hosts_they_infect):
        self.records = records_obj.records
        self.types = types
        self.viruses_and_hosts_they_infect = viruses_and_hosts_they_infect
        self.main_attributes_all_species = records_obj.df_main_attributes_for_all_species


    # get dictionary of tuples: {{'type1':(mean1,std1)},{'type2':(mean2,std2)}....}
    # helpful link: https://pythonforundergradengineers.com/python-matplotlib-error-bars.html
    def bar_chart_histogram(self,dict):
        #self.records = None
        # get list with size by type (parameter)
        # send the lists to the histogram with color (parameter) and column name (parameter)
        types = dict.keys()
        x_pos = np.arange(len(types)) # number of bars
        means = [x[0] for x in dict.values()]  # means or heights of the bars
        errors = [x[1] for x in dict.values()] # the error bars and the standard deviations

        # Build the plot
        fig, ax = plt.subplots()
        ax.bar(x_pos, means, yerr=errors, align='center', alpha=0.5, ecolor='black', capsize=10,color=['blue', 'red', 'white', 'white'],edgecolor="black")
        ax.set_ylabel('Genome size')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(types)
        ax.set_title('(B) Average genome sizes of Podoviruses, Myoviruses, Prochlorococcus and Synechococcus')
        ax.yaxis.grid(True)

        # Save the figure and show
        plt.tight_layout()
        plt.savefig('figure1_B_bar_plot_with_error_bars.png')
        plt.show()

    def get_mean_and_std(self):
        mean_and_std_of_types = {}
        for organism_type in self.types:
            mean_and_std_of_types[organism_type] = self.get_mean_and_std_by_type(organism_type)
        return mean_and_std_of_types

    def get_mean_and_std_by_type(self, organism_type):
        sizes = numpy.array([])
        for record in self.records:
            if organism_type in record.taxonomy:
                sizes = numpy.append(sizes, record.main_attributes_dictionary['genome_size'])
        return sizes.mean(), sizes.std()

    #-----

    def scatter_plot(self):
        # Synechococcus organisms:
        syne_gc_percentage = self.main_attributes_all_species.loc[self.main_attributes_all_species['family'] == 'Synechococcus', '%_GC_in_genes'].tolist()
        syne_trna = self.main_attributes_all_species.loc[self.main_attributes_all_species['family'] == 'Synechococcus', 'tRNA'].tolist()
        # LL-Prochlorococcus organisms:
        llproch_gc_percentage = self.main_attributes_all_species.loc[self.main_attributes_all_species['family'] == 'LL-Prochlorococcus', '%_GC_in_genes'].tolist()
        llproch_trna = self.main_attributes_all_species.loc[self.main_attributes_all_species['family'] == 'LL-Prochlorococcus', 'tRNA'].tolist()
        # HL-Prochlorococcus organisms:
        hlproch_gc_percentage = self.main_attributes_all_species.loc[self.main_attributes_all_species['family'] == 'HL-Prochlorococcus', '%_GC_in_genes'].tolist()
        hlproch_trna = self.main_attributes_all_species.loc[self.main_attributes_all_species['family'] == 'HL-Prochlorococcus', 'tRNA'].tolist()
        # Myoviridae organisms:
        myo_df = self.main_attributes_all_species.loc[self.main_attributes_all_species['family'] == 'Myoviridae']
        # Myoviridae (host: Synechococcus):
        myo_infect_syne_gc_percentage=[]
        myo_infect_syne_trna=[]
        # Myoviridae (host: Prochlorococcus):
        myo_infect_llproch_gc_percentage=[]
        myo_infect_llproch_trna=[]
        myo_infect_hlproch_gc_percentage=[]
        myo_infect_hlproch_trna=[]

        for index, row in myo_df.iterrows():
            myo_host_list = self.viruses_and_hosts_they_infect.get(row['name'])
            if any("Synechococcus" in s for s in myo_host_list):
                myo_infect_syne_gc_percentage.append(row['%_GC_in_genes'])
                if not isinstance(row['tRNA'], int):
                    row['tRNA']=0
                myo_infect_syne_trna.append(row['tRNA'])
            if any("LL-Prochlorococcus"in s for s in myo_host_list):
                myo_infect_llproch_gc_percentage.append(row['%_GC_in_genes'])
                if not isinstance(row['tRNA'], int):
                    row['tRNA']=0
                myo_infect_llproch_trna.append(row['tRNA'])
            if any("HL-Prochlorococcus"in s for s in myo_host_list):
                myo_infect_hlproch_gc_percentage.append(row['%_GC_in_genes'])
                if not isinstance(row['tRNA'], int):
                    row['tRNA']=0
                myo_infect_hlproch_trna.append(row['tRNA'])

        # Podoviridae organisms:
        podo_df = self.main_attributes_all_species.loc[self.main_attributes_all_species['family'] == 'Podoviridae']
        # Podoviridae (host: Synechococcus):
        podo_infect_syne_gc_percentage = []
        podo_infect_syne_trna = []
        # Podoviridae (host: Prochlorococcus):
        podo_infect_llproch_gc_percentage = []
        podo_infect_llproch_trna = []
        podo_infect_hlproch_gc_percentage = []
        podo_infect_hlproch_trna = []

        for index, row in podo_df.iterrows():
            podo_host_list = self.viruses_and_hosts_they_infect.get(row['name'])
            if any("Synechococcus" in s for s in podo_host_list):
                podo_infect_syne_gc_percentage.append(row['%_GC_in_genes'])
                if not isinstance(row['tRNA'], int):
                    row['tRNA']=0
                podo_infect_syne_trna.append(row['tRNA'])
            if any("LL-Prochlorococcus" in s for s in podo_host_list):
                podo_infect_llproch_gc_percentage.append(row['%_GC_in_genes'])
                if not isinstance(row['tRNA'], int):
                    row['tRNA']=0
                podo_infect_llproch_trna.append(row['tRNA'])
            if any("HL-Prochlorococcus" in s for s in podo_host_list):
                podo_infect_hlproch_gc_percentage.append(row['%_GC_in_genes'])
                if not isinstance(row['tRNA'], int):
                    row['tRNA']=0
                podo_infect_hlproch_trna.append(row['tRNA'])

        plt.scatter(syne_gc_percentage, syne_trna, c="none",marker="o",edgecolor="black")

        plt.scatter(llproch_gc_percentage, llproch_trna, c="none",marker="^",edgecolor="black")
        plt.scatter(hlproch_gc_percentage, hlproch_trna, c="none",marker="^",edgecolor="black")

        plt.scatter(myo_infect_syne_gc_percentage, myo_infect_syne_trna, c="red",linewidths=2,marker="o",edgecolor="none")
        plt.scatter(myo_infect_llproch_gc_percentage, myo_infect_llproch_trna, c="red",linewidths=2,marker="^",edgecolor="none")
        plt.scatter(myo_infect_hlproch_gc_percentage, myo_infect_hlproch_trna, c="red",linewidths=2,marker="^",edgecolor="none")

        plt.scatter(podo_infect_syne_gc_percentage, podo_infect_syne_trna, c="b",marker="o",edgecolor="none")
        plt.scatter(podo_infect_llproch_gc_percentage, podo_infect_llproch_trna, c="b",marker="^",edgecolor="none")
        plt.scatter(podo_infect_hlproch_gc_percentage, podo_infect_hlproch_trna, c="b",marker="^",edgecolor="none")

        plt.title("(A) The tRNA Gene Copy Number of various cyanobacteria and cyanophages against their GC Content")
        plt.legend(["Synechococcus", "LL-Prochlorococcus","HL-Prochlorococcus", "Myo (host: Synechococcus)", "Myo (host: LL-Prochlorococcus)", "Myo (host: HL-Prochlorococcus)","Podo (host: Synechococcus)", "Podo (host: LL-Prochlorococcus)", "Podo (host: HL-Prochlorococcus)"])

        plt.xlabel("GC content (%)")
        plt.ylabel("tRNA Gene Copy Number")
        plt.show()