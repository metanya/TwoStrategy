import sys

from Bio.Seq import reverse_complement
from Bio.SeqRecord import SeqRecord


# Copyright(C) 2009 Iddo Friedberg & Ian MC Fleming
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment


def get_interregions(genbank_record, intergene_length=1):
    size = 0
    gene_list_plus = []
    gene_list_minus = []
    intergenic_records = []

    get_gene_features(genbank_record, gene_list_minus, gene_list_plus)
    if len(gene_list_plus) != 0:
        size += get_intergene(genbank_record, gene_list_plus, intergene_length, intergenic_records, size, "+")
    else:
        intergenic_records.append(
            SeqRecord(
                genbank_record.seq,
                id="%s-ign-%d" % (genbank_record.name, 0),
                description="%s %d-%d %s"
                            % (genbank_record.name, len(genbank_record.seq), 0, "+"),
            ))
    if len(gene_list_minus) != 0:
        size = get_intergene(genbank_record, gene_list_minus, intergene_length, intergenic_records, size, "-")
    else:
        intergenic_records.append(
            SeqRecord(
                reverse_complement(genbank_record.seq),
                id="%s-ign-%d" % (genbank_record.name, 0),
                description="%s %d-%d %s"
                            % (genbank_record.name, len(genbank_record.seq), 0, "-"),
            ))
        size += len(genbank_record.seq)
    print(size)
    return intergenic_records, size


def get_gene_features(genbank_record, gene_list_minus, gene_list_plus):
    for feature in genbank_record.features:
        if feature.type == "gene" and feature.location_operator is None:
            mystart = feature.location.start.position
            myend = feature.location.end.position
            if feature.strand == -1:
                gene_list_minus.append((mystart, myend, -1))
            elif feature.strand == 1:
                gene_list_plus.append((mystart, myend, 1))
            else:
                sys.stderr.write("No strand indicated %d-%d. Assuming +\n" % (mystart, myend))
                gene_list_plus.append((mystart, myend, 1))


def get_intergene(genbank_record, gene_list, intergene_length, intergenic_records, size, strand):
    for i, pospair in enumerate(gene_list):
        if i - 1 < 0:
            last_end = 0
        else:
            last_end = gene_list[i - 1][1]
        this_start = pospair[0]
        if this_start - last_end >= intergene_length:
            add_intergenic(genbank_record, i, intergenic_records, last_end, strand, this_start)
        size += this_start - last_end
    if len(gene_list) != 0:
        add_intergenic(genbank_record, len(gene_list), intergenic_records, len(genbank_record.seq), strand,
                       gene_list[len(gene_list) - 1][1])
        size += len(genbank_record.seq) - gene_list[len(gene_list) - 1][1]
    return size


def add_intergenic(genbank_record, i, intergenic_records, last_end, strand, this_start):
    intergene_seq = genbank_record.seq[last_end:this_start]
    if strand == "-":
        intergene_seq = reverse_complement(intergene_seq)
    strand_string = strand
    intergenic_records.append(
        SeqRecord(
            intergene_seq,
            id="%s-ign-%d" % (genbank_record.name, i),
            description="%s %d-%d %s"
                        % (genbank_record.name, last_end + 1, this_start, strand_string),
        )
    )


def get_protein_gc_number(feature, seq):
    seq = seq[feature.location.start.position: feature.location.end.position]
    counter_g = seq.count('G')
    counter_c = seq.count('C')
    return counter_c + counter_g


def get_protein_gc_percentage(feature, seq):
    seq = seq[feature.location.start.position: feature.location.end.position]
    counter_g = seq.count('G')
    counter_c = seq.count('C')
    gc_percentage = ((counter_c + counter_g) / len(seq)) * 100
    return gc_percentage


def positive_negative_amino_acids(translation):
    trans_length = len(translation)
    counter_e = translation.count('E')
    counter_d = translation.count('D')
    counter_k = translation.count('K')
    counter_h = translation.count('H')
    counter_r = translation.count('R')
    positive = ((counter_r + counter_h + counter_k) / trans_length) * 100
    negative = ((counter_e + counter_d) / trans_length) * 100
    return positive, negative


def create_data_dictionary(record_gb):
    features = {}
    seq = record_gb.seq.upper()
    # record_gb.get_interregions()
    for feature in record_gb.features:
        trans_table = ""
        translation = ""
        gene = ""
        locus_tag = ""
        positive = -1
        negative = -1
        product = ""
        if 'transl_table' in feature.qualifiers:
            trans_table = feature.qualifiers['transl_table'][0]
        if 'translation' in feature.qualifiers:
            translation = feature.qualifiers['translation'][0]
            positive, negative = positive_negative_amino_acids(translation)
        if 'gene' in feature.qualifiers:
            gene = feature.qualifiers['gene'][0]
        if 'locus_tag' in feature.qualifiers:
            locus_tag = feature.qualifiers['locus_tag'][0]
        if 'product' in feature.qualifiers:
            product = feature.qualifiers['product'][0]

        gc_number = get_protein_gc_number(feature, seq)
        gc_percentage = get_protein_gc_percentage(feature, seq)
        name = str(feature.location.start) + " " + feature.type
        features[name] = {"start": feature.location.start,
                          "end": feature.location.end,
                          "length": feature.location.end - feature.location.start,
                          "type": feature.type,
                          "trans_table": trans_table,
                          "translation": translation,
                          "gene": gene,
                          "locus_tag": locus_tag,
                          "gc_number": gc_number,
                          "gc_percentage": gc_percentage,
                          "strand": feature.location.strand,
                          "amino_acid_positive_percentage": positive,
                          "amino_acid_negative_percentage": negative,
                          "product": product}

    return features
