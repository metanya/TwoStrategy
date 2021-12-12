import matplotlib.pyplot as plt
import numpy
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter


class Figures:
    def __init__(self, records, types):
        self.records = records
        self.types = types

    def histogram(self):
        self.records = None
        # get list with size by type (parameter)
        # send the lists to the histogram with color (parameter) and column name (parameter)

    def get_mean(self):
        mean_of_types = {}
        for organism_type in self.types:
            mean_of_types[organism_type] = self.get_mean_by_type(organism_type)
        return mean_of_types

    def get_mean_by_type(self, organism_type):
        sizes = numpy.array([])
        for record in self.records:
            if organism_type in record.taxonomy:
                sizes = numpy.append(sizes, record.main_attributes_dictionary['genome_size'])
        return sizes.mean(), sizes.std()
