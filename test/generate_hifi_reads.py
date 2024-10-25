import sys
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def introduce_errors(sequence, error_rate):
    num_errors = int(len(sequence) * error_rate)
    error_positions = random.sample(range(len(sequence)), num_errors)
    sequence_list = list(sequence)
    for pos in error_positions:
        sequence_list[pos] = random.choice('ATCG')
    return ''.join(sequence_list)

def generate_hifi_reads(input_fasta, output_fasta, coverage, read_length, error_rate):
    genome_length = 0
    sequences = []

    # Read the input FASTA file
    for record in SeqIO.parse(input_fasta, "fasta"):
        genome_length += len(record.seq)
        sequences.append(str(record.seq))

    num_reads = int((genome_length * coverage) / read_length)

    with open(output_fasta, 'w') as output_handle:
        for seq in sequences:
            seq_length = len(seq)
            num_reads_for_seq = int((seq_length * coverage) / read_length)

            for i in range(num_reads_for_seq):
                start = random.randint(0, seq_length - read_length)
                read_sequence = seq[start:start + read_length]
                read_sequence_with_errors = introduce_errors(read_sequence, error_rate)

                # Randomly decide if the read should be reverse-complemented
                if random.choice([True, False]):
                    read_sequence_with_errors = str(Seq(read_sequence_with_errors).reverse_complement())

                read_id = f"read_{i+1}"
                read_record = SeqRecord(Seq(read_sequence_with_errors), id=read_id, description="")
                SeqIO.write(read_record, output_handle, "fasta-2line")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_hifi_reads.py <input_fasta> <output_fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    coverage = 50
    read_length = 10000
    error_rate = 0.0001  # 0.01% error rate

    generate_hifi_reads(input_fasta, output_fasta, coverage, read_length, error_rate)

