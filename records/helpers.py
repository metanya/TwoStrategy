def create_main_attributes_dictionary_for_all_species(records):
    all_keys = sorted(set().union(*(d.main_attributes_dictionary.keys() for d in records)))
    main_attributes_dictionary_for_all_species = {}
    for record in records:
        new_dictionary = get_attributes(all_keys, record)
        main_attributes_dictionary_for_all_species[record.record_id] = new_dictionary
    return main_attributes_dictionary_for_all_species


def get_attributes(all_keys, record):
    new_dictionary = {}
    for attribute in all_keys:
        if attribute in record.main_attributes_dictionary.keys():
            new_dictionary[attribute] = record.main_attributes_dictionary[attribute]
        else:
            new_dictionary[attribute] = ""
    return new_dictionary
