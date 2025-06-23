
from Bio import SeqIO, Align

def read_fasta(path):
    try:
        record=SeqIO.read(path,"fasta")
        return str(record.seq)
    except Exception as e:
        print(f"Error opening path: {e}")
        exit()
    


def compute_alignment_stats(alignment):
    formatted_alignment=alignment.format().split('\n')
    aligned_seq1=formatted_alignment[0].strip()
    aligned_seq2=formatted_alignment[2].strip()

    matches = 0
    gaps = 0
    aligned_length = 0

    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == '-' or b == '-':
            gaps += 1
        else:
            aligned_length += 1
            if a == b:
                matches += 1

    identity_pct = (matches / aligned_length) * 100 if aligned_length > 0 else 0
    gap_pct = (gaps / len(aligned_seq1)) * 100

    return {
        "identity (%)": round(identity_pct, 2),
        "gap (%)": round(gap_pct, 2),
        "aligned_length": aligned_length
    }

def alignment(seq1,seq2,mode):
    try:
        aligner=Align.PairwiseAligner()
        aligner.mode=mode
        aligner.match_score=2
        aligner.mismatch_score=-1
        aligner.open_gap_score=-2
        aligner.extend_gap_score=-0.5

        alignments=aligner.align(seq1,seq2)
        top_alignment=alignments[0]

        return top_alignment
    except Exception as e:
        print(f"Excepion: {e}")
        exit()



path1=input("Enter the fasta path of your query:").strip()
path2=input("Enter the fasta path of your target:").strip()



seq1=read_fasta(path1)
seq2=read_fasta(path2)



mode= (input("Enter the mode of pairwise alignment, 'global' or 'local':")).strip().lower()
if mode not in ['global','local']:
    print ("Inavalid mode")
    exit()


top_alignments= alignment(seq1,seq2,mode)
score= top_alignments.score
stats=compute_alignment_stats(top_alignments)

print(f"The top alignment for the {mode.capitalize()} alignment is:")
print(top_alignments)
print(f"The score for top the alignment is: {score}")
print(f"The percentage of identity is {stats['identity (%)']}%")
print(f"The percentage of gap is {stats['gap (%)']}%")
print(f"The aligned length is {stats['aligned_length']}")


