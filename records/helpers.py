# data frame for all records together
from future.types import newint


def create_main_attributes_dictionary_for_all_species(records):
    all_keys = set().union(*(d.main_attributes_dictionary.keys() for d in records))
    main_attributes_dictionary_for_all_species = {}

    for record in records:
        new_dictionary = {}
        for attribute in all_keys:
            if attribute in record.main_attributes_dictionary.keys():
                new_dictionary[attribute] = record.main_attributes_dictionary[attribute]
            else:
                new_dictionary[attribute] = ""
        main_attributes_dictionary_for_all_species[record.record_id] = new_dictionary
        # main_attributes_dictionary_for_all_species[record.record_id] = {
        #     "gene_num": record.main_attributes_dictionary["gene"],
        #     "cds_num": record.main_attributes_dictionary["CDS"],
        #     "trna_num": record.main_attributes_dictionary["tRNA"],
        #     # "genome_size": record.main_attributes_dictionary["genome_size"],
        #     # "percentage_of_genes_in_genome": record.main_attributes_dictionary["percentage_of_genes_in_genome"],
        #     # "percentage_of_intergene_in_genome": record.main_attributes_dictionary["percentage_of_intergene_in_genome"],
        #     # "percentage_of_GC_in_genome": record.main_attributes_dictionary["percentage_of_GC_in_genome"],
        #     # "percentage_of_GC_in_genes": record.main_attributes_dictionary["percentage_of_GC_in_genes"],
        #     # "percentage_of_GC_in_intergene": record.main_attributes_dictionary["percentage_of_GC_in_intergene"]
        # }
    print(main_attributes_dictionary_for_all_species)
    return main_attributes_dictionary_for_all_species
