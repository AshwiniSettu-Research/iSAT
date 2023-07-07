from flask import Flask, render_template, request
import io
import base64
import json
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from colorama import init, Fore, Style
from termcolor import colored
from scipy.sparse import lil_matrix

app = Flask(__name__)


def calculate_fingerprint(protein_sequence):
    """
    Calculates the structure-based fingerprint of a protein sequence.

    Args:
        protein_sequence: The protein sequence.

    Returns:
        The structure-based fingerprint of the protein sequence.
    """

    # Initialize the fingerprint.
    fingerprint = np.zeros(20)

    # Loop over the amino acids in the protein sequence.
    for amino_acid in protein_sequence:
        # Calculate the weight of the amino acid.
        weight = get_weight(amino_acid)

        # Update the fingerprint.
        fingerprint += weight

    # Return the fingerprint.
    return fingerprint


def get_weight(amino_acid):
    """
    Gets the weight of an amino acid.

    Args:
        amino_acid: The amino acid.

    Returns:
        The weight of the amino acid.
    """

    # Get the weight of the amino acid from a dictionary.
    weights = {
        'A': 1,
        'C': 2,
        'D': 3,
        'E': 4,
        'F': 5,
        'G': 6,
        'H': 7,
        'I': 8,
        'L': 9,
        'K': 10,
        'M': 11,
        'N': 12,
        'O': 13,
        'P': 14,
        'Q': 15,
        'R': 16,
        'S': 17,
        'T': 18,
        'U': 19,
        'V': 20,
        'W': 21,
        'X': 22,
        'Y': 23
    }

    # Return the weight of the amino acid.
    return weights[amino_acid]


def sfa(a, b, w_m, w_d, g, alignment_matrix=None):
    """
    Aligns two protein sequences using the SFA algorithm.

    Args:
        a: The first protein sequence.
        b: The second protein sequence.
        w_m: The weight of a match.
        w_d: The weight of a mismatch.
        g: The gap penalty.
        alignment_matrix: The alignment matrix (optional).

    Returns:
        The similarity score, accuracy, alignment, and traceback alignment.
    """

    # Check if the protein sequences contain any decimal characters.
    if any(c.isdecimal() for c in a):
        raise ValueError('The protein sequence cannot contain any decimal characters.')
    if any(c.isdecimal() for c in b):
        raise ValueError('The protein sequence cannot contain any decimal characters.')

    # Calculate the structure-based fingerprints of the proteins.
    fp_a = calculate_fingerprint(a)
    fp_b = calculate_fingerprint(b)

    # Initialize the alignment matrix if not provided.
    if alignment_matrix is None:
        alignment_matrix = np.zeros((len(a) + 1, len(b) + 1))

    # Fill in the alignment matrix.
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                alignment_matrix[i][j] = alignment_matrix[i - 1][j - 1] + w_m
            else:
                alignment_matrix[i][j] = max(alignment_matrix[i - 1][j - 1] + w_d,
                                             alignment_matrix[i - 1][j] + g,
                                             alignment_matrix[i][j - 1] + g)

    # Traceback through the alignment matrix to identify the alignment.
    alignment = []
    i = len(a)
    j = len(b)
    while i > 0 or j > 0:
        if a[i - 1] == b[j - 1]:
            alignment.append((a[i - 1], 'match'))
            i -= 1
            j -= 1
        else:
            if alignment_matrix[i - 1][j - 1] + w_d >= max(alignment_matrix[i - 1][j] + g,
                                                          alignment_matrix[i][j - 1] + g):
                alignment.append((a[i - 1], 'mismatch'))
                i -= 1
                j -= 1
            else:
                if alignment_matrix[i - 1][j] + g > alignment_matrix[i][j - 1] + g:
                    alignment.append((a[i - 1], 'deletion'))
                    i -= 1
                else:
                    alignment.append((b[j - 1], 'insertion'))
                    j -= 1

    # Calculate the similarity score.
    similarity_score = alignment_matrix[len(a)][len(b)]

    # Calculate the accuracy.
    accuracy = similarity_score / len(a)

    # Create a traceback alignment with colors.
    traceback_alignment = []
    for amino_acid, operation in alignment:
        if operation == 'match':
            traceback_alignment.append((amino_acid, 'black'))
        elif operation == 'mismatch':
            traceback_alignment.append((amino_acid, 'red'))
        elif operation == 'deletion':
            traceback_alignment.append((amino_acid, 'blue'))
        elif operation == 'insertion':
            traceback_alignment.append((amino_acid, 'green'))

    # Print the alignment in color.
    print('Alignment:')
    for amino_acid, operation in alignment:
        print(f'{amino_acid:>4} {operation:>4}')

    # Create a heatmap of the alignment matrix.
    heatmap = plt.imshow(alignment_matrix, cmap='hot')
    plt.colorbar()
    plt.title('Alignment Matrix Heatmap')
    plt.show()

    # Create a scatter plot of the aligned sequences.
    aligned_a = [amino_acid for amino_acid, _ in alignment]
    aligned_b = [amino_acid for _, amino_acid in alignment]
    plt.scatter(range(len(aligned_a)), aligned_a, c=[color for _, color in traceback_alignment])
    plt.scatter(range(len(aligned_b)), aligned_b, c=[color for _, color in traceback_alignment])
    plt.xlabel('Alignment Position')
    plt.ylabel('Aligned Amino Acids')
    plt.show()

    # Print the traceback alignment as a single line.
    traceback_alignment_seq = ''.join([amino_acid for amino_acid, _ in traceback_alignment[::-1]])
    print(f"Traceback alignment: {traceback_alignment_seq}")

    return similarity_score, accuracy, alignment, traceback_alignment, traceback_alignment_seq, alignment_matrix

def sequence_consensus(aligned_sequences):
    consensus = ''
    for positions in zip(*aligned_sequences):
        counts = {}
        for amino_acid in positions:
            if amino_acid in counts:
                counts[amino_acid] += 1
            else:
                counts[amino_acid] = 1
        consensus += max(counts, key=counts.get)
    return consensus

# Define the route for the home page
@app.route('/')
def index():
    return render_template('index.html')

# Define the route to process the form submission
@app.route('/process', methods=['POST'])
def process_form():
    # Retrieve the input values from the form
    protein_sequence_a = request.form['protein_sequence_a']
    protein_sequence_b = request.form['protein_sequence_b']
    

    # Call your existing function or scripts here to process the input and obtain the output
    similarity_score, accuracy, alignment, traceback_alignment, traceback_alignment_seq, alignment_matrix = sfa(protein_sequence_a, protein_sequence_b, 1, -1, -2)

    # Create a table-like output as a string
    output_table = 'Alignment Results:\n'
    output_table += '---------------------------------------------------------------------------------------------------------------------------------------------------------------------\n'
    output_table += f'Length of Sequence 1: {len(protein_sequence_a)}\n'
    output_table += 'Sequence 1:\n\n'
    for i in range(0, len(protein_sequence_a), 160):  # Print 80 characters per row
        output_table += protein_sequence_a[i:i+160] + '\n'
    output_table += f'Length of Sequence 2: {len(protein_sequence_b)}\n'
    output_table += 'Sequence 2:\n\n'
    for i in range(0, len(protein_sequence_b), 160):  # Print 80 characters per row
        output_table += protein_sequence_b[i:i+160] + '\n'
    output_table += '---------------------------------------------------------------------------------------------------------------------------------------------------------------------\n'
    output_table += f'Traceback alignment: {traceback_alignment_seq}\n\n'   
    output_table += '---------------------------------------------------------------------------------------------------------------------------------------------------------------------\n'
    output_table += f'Similarity score: {similarity_score}\n'
    output_table += f'Accuracy: {accuracy:.2%}\n\n'
    output_table += '----------------------------------------------------------------------------------------------------------------------------------------------------------------------\n'

  # Create a heatmap of the alignment matrix.
    plt.figure()
    plt.imshow(alignment_matrix, cmap='hot')
    plt.colorbar()
    plt.title('Alignment Matrix Heatmap')
    heatmap_image = io.BytesIO()
    plt.savefig(heatmap_image, format='png')
    heatmap_image.seek(0)
    heatmap_base64 = base64.b64encode(heatmap_image.getvalue()).decode('utf-8')

    # Create a scatter plot of the aligned sequences.
    plt.figure()
    aligned_a = [amino_acid for amino_acid, _ in alignment]
    aligned_b = [amino_acid for _, amino_acid in alignment]
    plt.scatter(range(len(aligned_a)), aligned_a, c=[color for _, color in traceback_alignment])
    plt.scatter(range(len(aligned_b)), aligned_b, c=[color for _, color in traceback_alignment])
    plt.xlabel('Alignment Position')
    plt.ylabel('Aligned Amino Acids')
    scatter_plot_image = io.BytesIO()
    plt.savefig(scatter_plot_image, format='png')
    scatter_plot_image.seek(0)
    scatter_plot_base64 = base64.b64encode(scatter_plot_image.getvalue()).decode('utf-8')

    # Create the alignment table.
    alignment_table = [('Position','Sequence','Correct','Mismatches','Insertions','Deletions')]
    for i, (amino_acid, operation) in enumerate(alignment, start=1):
        position_number = i
        sequence_letter = amino_acid
        correct_letter = amino_acid if operation == 'match' else '-'
        mismatches = amino_acid if operation == 'mismatch' else '-'
        insertions = amino_acid if operation == 'insertion' else '-'
        deletions = amino_acid if operation == 'deletion' else '-'
        alignment_table.append((position_number, sequence_letter, correct_letter, mismatches, insertions, deletions))

    consensus_sequence = sequence_consensus([aligned_a, aligned_b])

    return render_template('result.html', output_table=output_table, alignment_table=alignment_table,
                       heatmap_base64=heatmap_base64, scatter_plot_base64=scatter_plot_base64,
                       consensus_sequence=consensus_sequence)

if __name__ == '__main__':
    app.run(debug=True)
