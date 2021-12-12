import numpy
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import colors
from matplotlib.ticker import PercentFormatter


class Figures:
    def __init__(self, records, types):
        self.records = records
        self.types = types

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
        ax.bar(x_pos, means, yerr=errors, align='center', alpha=0.5, ecolor='black', capsize=10)
        ax.set_ylabel('Genome size (pb)')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(types)
        ax.set_title('(B) Average genome sizes of podoviruses,myoviruses, Prochlorococcus and Synechococcus')
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
