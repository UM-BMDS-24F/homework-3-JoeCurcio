from Bio import SeqIO
from Bio.Blast import NCBIXML
import subprocess
import os
import io  

# File paths
human_file = "human.fa"
mouse_file = "mouse.fa"
mouse_db = "mouse_db"  

# Parsing human sequences
human_sequences = list(SeqIO.parse(human_file, "fasta"))

# Creating a BLAST database 
if not os.path.exists(mouse_db + ".pin"):
    subprocess.run(["makeblastdb", "-in", mouse_file, "-dbtype", "prot", "-out", mouse_db])

# Running BLAST for each human sequence
output_file = "blast_results.txt"
with open(output_file, "w") as f:
    for human_seq in human_sequences:
        with open("temp_human_query.fa", "w") as temp_query_file:
            SeqIO.write(human_seq, temp_query_file, "fasta")
        blast_command = [
            "blastp", 
            "-query", "temp_human_query.fa", 
            "-db", mouse_db, 
            "-evalue", "1e-5", 
            "-outfmt", "5", 
            "-num_alignments", "1"
        ]
        
        # Executing BLAST
        result = subprocess.run(blast_command, capture_output=True, text=True)
        
        # Checking for errors
        if result.returncode != 0:
            print(f"Error running BLAST: {result.stderr}")
            continue

        blast_record = NCBIXML.read(io.StringIO(result.stdout))

        # Processing results if alignments are found
        if blast_record.alignments:
            top_alignment = blast_record.alignments[0]
            top_hsp = top_alignment.hsps[0]
            
            # Extracting the Mouse Sequence ID from the hit ID
            mouse_id = top_alignment.hit_id.split("|")[1]
            
            # Writing the results to the output file
            f.write(f"Human Seq ID: {human_seq.id}\n")
            f.write(f"Mouse Seq ID: {mouse_id}\n")
            f.write(f"Alignment:\n{top_hsp.sbjct}\n")
            f.write(f"E-value: {top_hsp.expect}\n")
            f.write(f"Bitscore: {top_hsp.bits}\n\n")
