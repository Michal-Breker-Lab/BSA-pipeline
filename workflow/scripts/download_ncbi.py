from snakemake.script import snakemake
from pathlib import Path
import subprocess
import tempfile
import zipfile
import shutil
from Bio import Entrez
import subprocess
import time
import sys


log = open(snakemake.log[0], "w")
sys.stderr = sys.stdout = log


def download_assembly(assembly_acc, out_dir):
    """Download genome or annotation from NCBI using datasets command line tool."""
    out_dir = Path(out_dir)
    if not assembly_acc.startswith(("GCF_", "GCA_")):
        raise ValueError("Assembly accession must start with GCF_ or GCA_")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        subprocess.run(
            [
                "datasets",
                "download",
                "genome",
                "accession",
                assembly_acc,
                "--filename",
                f"{tmpdir}/assembly.zip",
                "--no-progressbar",
                "--include",
                "genome,gff3"
            ],
            check=True,
        )

        with zipfile.ZipFile(f"{tmpdir}/assembly.zip", 'r') as zip_ref:
            zip_ref.extractall(tmpdir)

        assembly_folder = Path(f"{tmpdir}/ncbi_dataset/data/{assembly_acc}/")

        genome_file = list(assembly_folder.glob(f"{assembly_acc}*.fna"))
        gff_file = list(assembly_folder.glob(f"*.gff"))

        if len(genome_file) != 1:
            raise FileNotFoundError(f"Expected one fasta file for {assembly_acc}, found {len(genome_file)}")
        shutil.move(genome_file[0], out_dir / f"{assembly_acc}.fa")

        gff_out_path = out_dir / f"{assembly_acc}.gff"
        if len(gff_file) == 1:
            shutil.move(gff_file[0], gff_out_path)
        else:
            print("Warning:{assembly_acc} GFF file not found")
            Path(gff_out_path).touch()


def get_content(acc, rettype):
    fetch_handle = Entrez.efetch(db="nuccore", id=acc, rettype=rettype, retmode="text")
    content = fetch_handle.read()
    fetch_handle.close()
    if not content.strip():
        raise RuntimeError(f"Entrez Efetch returned an empty response for ID {acc}.")
    
    return content


def download_sequence(accession_id, out_dir):
    out_dir = Path(out_dir)
    fasta_content = get_content(accession_id, "fasta")
    if fasta_content.startswith(">"):
        with open(out_dir / f"{accession_id}.fa", "w") as f:
            f.write(fasta_content.strip() + "\n")
    else:
        raise ValueError(f"Fasta file for {accession_id} is not valid")
    time.sleep(1)

    gff_content = get_content(accession_id, "gff")
    with open(out_dir / f"{accession_id}.gff", "w") as f:
        if gff_content.startswith("##gff-version"):
            f.write(gff_content.strip() + "\n")
        else:
            print(f"Warning: {accession_id} GFF file not found")
            f.write("")
    time.sleep(1)


Entrez.email = snakemake.params["email"]
out_dir = snakemake.params["out_dir"]

for a_acc in snakemake.params["ncbi_assemblies"]:
    download_assembly(a_acc, out_dir)
for s_acc in snakemake.params["ncbi_sequences"]:
    download_sequence(s_acc, out_dir)

