import sys

from Bio.Seq import reverse_complement
from Bio.SeqRecord import SeqRecord


# gene length = 178249


# Copyright(C) 2009 Iddo Friedberg & Ian MC Fleming
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment

def get_interregions(genbank_record, intergene_length=1):
    size = 0
    gene_list = []
    gene_list_plus = []
    gene_list_minus = []
    intergenic_records = []
    # Loop over the genome file, get the gene features on each of the strands
    for feature in genbank_record.features:
        if feature.type == "gene":
            mystart = feature.location.start.position
            myend = feature.location.end.position
            if feature.strand == -1:
                gene_list.append((mystart, myend, -1))
            elif feature.strand == 1:
                gene_list.append((mystart, myend, 1))
            else:
                sys.stderr.write(
                    "No strand indicated %d-%d. Assuming +\n" % (mystart, myend)
                )
                gene_list.append((mystart, myend, 1))
    for i, pospair in enumerate(gene_list[1:]):
        # Compare current start position to previous end position
        last_end = gene_list[i][1]
        this_start = pospair[0]
        if this_start - last_end >= intergene_length:
            intergene_seq = genbank_record.seq[last_end:this_start]
            if feature.strand == -1:
                strand_string = "-"
            elif feature.strand == 1:
                strand_string = "+"
            intergenic_records.append(
                SeqRecord(
                    intergene_seq,
                    id="%s-ign-%d" % (genbank_record.name, i),
                    description="%s %d-%d %s"
                                % (genbank_record.name, last_end + 1, this_start, strand_string),
                )
            )
            size += this_start - last_end
    print(size)
    return intergenic_records

    # print(size)
    # return intergenic_records
    # size = 0
    # gene_list_plus = []
    # gene_list_minus = []
    # intergenic_records = []
    # # Loop over the genome file, get the gene features on each of the strands
    # for feature in genbank_record.features:
    #     if feature.type == "gene":
    #         mystart = feature.location.start.position
    #         myend = feature.location.end.position
    #         if feature.strand == -1:
    #             gene_list_minus.append((mystart, myend, -1))
    #         elif feature.strand == 1:
    #             gene_list_plus.append((mystart, myend, 1))
    #         else:
    #             sys.stderr.write(
    #                 "No strand indicated %d-%d. Assuming +\n" % (mystart, myend)
    #             )
    #             gene_list_plus.append((mystart, myend, 1))
    #
    # for i, pospair in enumerate(gene_list_plus[1:]):
    #     # Compare current start position to previous end position
    #     last_end = gene_list_plus[i][1]
    #     this_start = pospair[0]
    #     if this_start - last_end >= intergene_length:
    #         intergene_seq = genbank_record.seq[last_end:this_start]
    #         strand_string = "+"
    #         intergenic_records.append(
    #             SeqRecord(
    #                 intergene_seq,
    #                 id="%s-ign-%d" % (genbank_record.name, i),
    #                 description="%s %d-%d %s"
    #                             % (genbank_record.name, last_end + 1, this_start, strand_string),
    #             )
    #         )
    #         size += this_start - last_end
    # for i, pospair in enumerate(gene_list_minus[1:]):
    #     last_end = gene_list_minus[i][1]
    #     this_start = pospair[0]
    #     if this_start - last_end >= intergene_length:
    #         intergene_seq = genbank_record.seq[last_end:this_start]
    #         intergene_seq = reverse_complement(intergene_seq)
    #         strand_string = "-"
    #         intergenic_records.append(
    #             SeqRecord(
    #                 intergene_seq,
    #                 id="%s-ign-%d" % (genbank_record.name, i),
    #                 description="%s %d-%d %s"
    #                             % (genbank_record.name, last_end + 1, this_start, strand_string),
    #             )
    #         )
    #         size += this_start - last_end
    # print(size)
    # return intergenic_records
