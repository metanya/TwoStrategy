from record.record import fields


def get_main_attributes(records):
    keys = sorted(fields)
    main_attributes = {}
    for record in records:
        main_attributes[record.record_id] = get_attributes(keys, record)
    return main_attributes


def get_attributes(keys, record):
    attributes = {}
    attributes_keys = record.main_attributes.keys()
    for key in keys:
        if key in attributes_keys:
            attributes[key] = record.main_attributes[key]
        else:
            attributes[key] = ""
    return attributes
