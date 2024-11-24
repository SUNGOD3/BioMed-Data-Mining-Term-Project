import numpy as np
from collections import defaultdict
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class RNAStructurePredictor:
    def __init__(self):
        # Define valid base pairs
        self.valid_pairs = {
            ('A', 'U'), ('U', 'A'),
            ('G', 'C'), ('C', 'G'),
            ('G', 'U'), ('U', 'G')
        }
        
    def read_fasta(self, filename):
        """Read sequences from FASTA file."""
        sequences = []
        with open(filename, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences.append(str(record.seq))
        return sequences

    def create_alignment(self, sequences):
        """Create a simple multiple sequence alignment."""
        # For this implementation, we assume sequences are pre-aligned
        # In a production environment, you would want to use a proper MSA tool
        max_len = max(len(seq) for seq in sequences)
        aligned_sequences = []
        
        for seq in sequences:
            if len(seq) < max_len:
                # Pad shorter sequences with gaps
                aligned_seq = seq + '-' * (max_len - len(seq))
            else:
                aligned_seq = seq
            aligned_sequences.append(aligned_seq)
            
        return aligned_sequences

    def calculate_covariation_matrix(self, alignment):
        """Calculate covariation scores between all positions."""
        length = len(alignment[0])
        covariation_matrix = np.zeros((length, length))
        
        for i in range(length):
            for j in range(i + 1, length):
                score = self._calculate_position_covariation(alignment, i, j)
                covariation_matrix[i, j] = score
                covariation_matrix[j, i] = score
                
        return covariation_matrix

    def _calculate_position_covariation(self, alignment, pos1, pos2):
        """Calculate covariation score between two positions."""
        pairs = defaultdict(int)
        valid_count = 0
        total_count = 0
        
        for seq in alignment:
            if seq[pos1] != '-' and seq[pos2] != '-':
                pairs[(seq[pos1], seq[pos2])] += 1
                total_count += 1
                if (seq[pos1], seq[pos2]) in self.valid_pairs:
                    valid_count += 1
                    
        if total_count == 0:
            return 0
            
        # Calculate covariation score based on valid base pairs frequency
        return valid_count / total_count if total_count > 0 else 0

    def predict_structure(self, alignment, min_score=0.75):
        """Predict secondary structure based on covariation analysis."""
        covariation_matrix = self.calculate_covariation_matrix(alignment)
        length = len(alignment[0])
        structure = ['.' for _ in range(length)]
        
        # Find base pairs with high covariation scores
        paired_positions = set()
        for i in range(length):
            if i in paired_positions:
                continue
                
            max_score = 0
            max_j = -1
            
            for j in range(i + 4, length):  # Minimum loop size of 3
                if j in paired_positions:
                    continue
                    
                if covariation_matrix[i, j] > max_score:
                    max_score = covariation_matrix[i, j]
                    max_j = j
                    
            if max_score >= min_score and max_j != -1:
                structure[i] = '('
                structure[max_j] = ')'
                paired_positions.add(i)
                paired_positions.add(max_j)
                
        return ''.join(structure)

    def format_output(self, sequences, structure):
        """Format the output as specified."""
        output = []
        for i, seq in enumerate(sequences):
            output.append(f">seq{i+1}")
            output.append(seq)
        output.append(structure)
        return '\n'.join(output)

def main(input_file):
    predictor = RNAStructurePredictor()
    
    # Read sequences
    sequences = predictor.read_fasta(input_file)
    
    # Create alignment
    alignment = predictor.create_alignment(sequences)
    
    # Predict structure
    structure = predictor.predict_structure(alignment)
    
    # Format and return output
    return predictor.format_output(alignment, structure)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python rna_structure_predictor.py <input_fasta_file>")
        sys.exit(1)
        
    result = main(sys.argv[1])
    print(result)
