# data frame for all records together
def create_main_attributes_dictionary_for_all_species(records):
    main_attributes_dictionary_for_all_species = {}
    for record in records:
        main_attributes_dictionary_for_all_species[record.record_id] = {
            "gene_num": record.main_attributes_dictionary["gene_num"],
            "cds_num": record.main_attributes_dictionary["cds_num"],
            "trna_num": record.main_attributes_dictionary["trna_num"],
            "genome_size": record.main_attributes_dictionary["genome_size"],
            "percentage_of_genes_in_genome": record.main_attributes_dictionary["percentage_of_genes_in_genome"],
            "percentage_of_intergene_in_genome": record.main_attributes_dictionary["percentage_of_intergene_in_genome"],
            "percentage_of_GC_in_genome": record.main_attributes_dictionary["percentage_of_GC_in_genome"],
            "percentage_of_GC_in_genes": record.main_attributes_dictionary["percentage_of_GC_in_genes"],
            "percentage_of_GC_in_intergene": record.main_attributes_dictionary["percentage_of_GC_in_intergene"]
        }
    print(main_attributes_dictionary_for_all_species)
    return main_attributes_dictionary_for_all_species
