import sys
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def calculate_protein_stats(sequence):
    """
    Calculate and return various protein sequence statistics.
    
    Args:
        sequence (str): Protein sequence (single-letter amino acid codes).
        
    Returns:
        dict: A dictionary of calculated protein statistics.
    """
    # Initialize the ProteinAnalysis object
    protein = ProteinAnalysis(sequence)

    # Amino acid composition
    amino_acid_composition = protein.count_amino_acids()

    # Molecular weight
    mol_weight = protein.molecular_weight()

    # Isoelectric point (pI)
    isoelectric_point = protein.isoelectric_point()

    # Extinction coefficient
    ext_coeff_reduced, ext_coeff_oxidized = protein.molar_extinction_coefficient()

    # Hydrophobicity
    hydrophobicity = protein.gravy()

    # Calculate Aliphatic Index
    aliphatic_index = calculate_aliphatic_index(sequence)

    # GRAVY (Grand Average of Hydropathy)
    gravy = protein.gravy()

    # Instability Index
    instability_index = protein.instability_index()

    # Flexibility
    flexibility = protein.flexibility()

    # Calculate Aromaticity
    aromaticity = protein.aromaticity()

    # Output all statistics as a dictionary
    protein_stats = {
        'Amino Acid Composition': amino_acid_composition,
        'Molecular Weight (Da)': mol_weight,
        'Isoelectric Point (pI)': isoelectric_point,
        'Extinction Coefficient (Reduced Cys)': ext_coeff_reduced,
        'Extinction Coefficient (Oxidized Cys)': ext_coeff_oxidized,
        'Hydrophobicity (GRAVY)': hydrophobicity,
        'Aliphatic Index': aliphatic_index,
        'Instability Index': instability_index,
        'Flexibility': flexibility,
        'Aromaticity': aromaticity
    }

    return protein_stats

def calculate_aliphatic_index(sequence):
    """
    Calculate the aliphatic index of a protein sequence.
    
    Args:
        sequence (str): Protein sequence (single-letter amino acid codes).
        
    Returns:
        float: The aliphatic index of the protein sequence.
    """
    # Aliphatic amino acids
    aliphatic_amino_acids = 'AVIL'
    total_residues = len(sequence)
    if total_residues == 0:
        return 0.0
    
    aliphatic_residues = sum(sequence.count(aa) for aa in aliphatic_amino_acids)
    return 100.0 * aliphatic_residues / total_residues

def print_protein_stats(header, stats):
    """
    Print protein statistics for a given header.
    
    Args:
        header (str): Header of the FASTA sequence.
        stats (dict): A dictionary of protein statistics.
    """
    print(f">{header}")
    for key, value in stats.items():
        print(f"{key}: {value}")
    print()  # Blank line for separation

def process_fasta_file(file_path):
    """
    Process a FASTA file to calculate and print statistics for each protein sequence.
    
    Args:
        file_path (str): Path to the FASTA file.
    """
    for record in SeqIO.parse(file_path, "fasta"):
        header = record.id
        sequence = str(record.seq)
        stats = calculate_protein_stats(sequence)
        print_protein_stats(header, stats)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 protein_analysis.py <fasta_file>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    process_fasta_file(fasta_file)
