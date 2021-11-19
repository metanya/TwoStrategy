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
