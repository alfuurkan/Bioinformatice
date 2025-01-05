import pandas as pd
from Bio.Align import PairwiseAligner

# Step 1: Load the dataset
def load_dataset(file_path):
    try:
        df = pd.read_csv(file_path, sep="\t")
        print("Dataset loaded successfully!")
        print("Columns in the dataset:", df.columns)  # Debug: Print column names
        return df
    except Exception as e:
        print(f"Error loading dataset: {e}")
        return None

# Step 2: Find the Best Matching GeneID
def find_best_matching_geneid(dataset, user_seq):
    if "GeneID" not in dataset.columns or "sequence" not in dataset.columns:
        print("Error: Required columns 'GeneID' and 'sequence' not found in the dataset.")
        return None

    aligner = PairwiseAligner()
    aligner.mode = 'global'  # Global alignment

    best_score = -float("inf")
    best_geneid = None
    best_alignment = None

    # Iterate through all sequences in the dataset
    for index, row in dataset.iterrows():
        target_sequence = row["sequence"]
        gene_id = row["GeneID"]
        alignment = aligner.align(user_seq, target_sequence)[0]
        score = alignment.score

        # Check if this alignment is the best so far
        if score > best_score:
            best_score = score
            best_geneid = gene_id
            best_alignment = alignment

    return best_geneid, best_score, best_alignment

# Step 3: Main Function
def main():
    # Load the dataset
    file_path = "covid.tsv"  # Update with your correct file path
    dataset = load_dataset(file_path)
    if dataset is None:
        return
    
    # User inputs a sequence
    user_sequence = input("Enter your sequence (A, T, C, G): ").strip().upper()
    if not user_sequence.isalpha() or any(base not in "ATCG" for base in user_sequence):
        print("Invalid sequence! Please use only A, T, C, and G.")
        return

    # Find the best matching GeneID
    print("\nFinding the best matching GeneID...")
    result = find_best_matching_geneid(dataset, user_sequence)
    if result is None:
        print("No matching sequence found.")
        return
    
    best_geneid, best_score, best_alignment = result

    # Display results
    print(f"\nBest Matching GeneID: {best_geneid}")
    print(f"Alignment Score: {best_score}")
    print("\nBest Alignment Result:")
    print(best_alignment)

if __name__ == "__main__":
    main()
