
import scrublet as scr
import pandas as pd
import scipy.io
import scipy.sparse
from pathlib import Path


def load_data(data_dir):
    """
    Load data (matrix, barcodes, features).
    
    Parameters:
    -----------
    data_dir : str
        Path to directory containing matrix.mtx, barcodes.tsv, features.tsv
    
    Returns:
    --------
    tuple
        (matrix, barcodes, features) - matrix is cells x genes for Scrublet
    """
    data_dir = Path(data_dir)
    print("Loading data...")
    
    # Load matrix (transpose to get cells x genes for Scrublet)
    matrix = scipy.io.mmread(data_dir / "matrix.mtx").T
    matrix = scipy.sparse.csr_matrix(matrix)
    
    # Load barcodes and features
    barcodes = pd.read_csv(data_dir / "barcodes.tsv", header=None)[0].values
    features = pd.read_csv(data_dir / "features.tsv", header=None)[0].values
    
    print(f"Loaded matrix: {matrix.shape} (cells x genes)")
    print(f"Barcodes: {len(barcodes)}")
    print(f"Features: {len(features)}")
    
    return matrix, barcodes, features


def run_scrublet_analysis(matrix, expected_doublet_rate=0.075):
    """
    Run Scrublet doublet detection.
    
    Parameters:
    -----------
    matrix : scipy.sparse matrix
        Expression matrix (cells x genes)
    expected_doublet_rate : float
        Expected doublet rate for the experiment
    
    Returns:
    --------
    tuple
        (scrub, doublet_scores, predicted_doublets)
    """
    print("Running Scrublet doublet detection...")
    
    scrub = scr.Scrublet(matrix, expected_doublet_rate=expected_doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    
    print(f"Doublet detection complete. Auto threshold: {scrub.threshold_:.4f}")
    return scrub, doublet_scores, predicted_doublets


def save_results(barcodes, doublet_scores, predicted_doublets, output_path):
    """
    Save doublet detection results to CSV.
    
    Parameters:
    -----------
    barcodes : np.ndarray
        Cell barcodes
    doublet_scores : np.ndarray
        Doublet scores
    predicted_doublets : np.ndarray
        Predicted doublet labels
    output_path : str
        Path to save the CSV file
    """
    # Ensure all arrays have the same length
    n_cells = len(doublet_scores)
    
    results_df = pd.DataFrame({
        'barcode': barcodes[:n_cells],
        'doublet_score': doublet_scores,
        'predicted_doublet': predicted_doublets
    })
    
    results_df.to_csv(output_path, index=False)
    print(f"Results saved to: {output_path}")


def export_filtered_data(matrix, barcodes, features, scrub, output_dir):
    """
    Export filtered data (doublets removed) 
    
    Parameters:
    -----------
    matrix : scipy.sparse matrix
        Original expression matrix (cells x genes)
    barcodes : np.ndarray
        Cell barcodes
    features : np.ndarray
        Gene names
    scrub : Scrublet object
        Fitted Scrublet object with called doublets
    output_dir : str
        Directory to save filtered data
    """
    # Filter out doublets
    clean_matrix = matrix[~scrub.predicted_doublets_]
    clean_barcodes = barcodes[~scrub.predicted_doublets_]
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Save filtered matrix (transpose back to genes x cells)
    scipy.io.mmwrite(output_dir / "matrix.mtx", clean_matrix.T)
    
    # Save filtered barcodes
    with open(output_dir / "barcodes.tsv", "w") as f:
        for barcode in clean_barcodes:
            f.write(f"{barcode}\n")
    
    # Copy original features (genes unchanged)
    with open(output_dir / "features.tsv", "w") as f:
        for feature in features:
            f.write(f"{feature}\n")
    
    print(f"Filtered data exported to: {output_dir}")
    return output_dir


def generate_final_report(matrix, scrub, expected_doublet_rate=0.075):
    """
    Generate a comprehensive analysis report.
    
    Parameters:
    -----------
    matrix : scipy.sparse matrix
        Original expression matrix (cells x genes)
    scrub : Scrublet object
        Fitted Scrublet object with called doublets
    expected_doublet_rate : float
        Expected doublet rate
    """
    total_cells = matrix.shape[0]
    doublets_removed = scrub.predicted_doublets_.sum()
    cells_remaining = total_cells - doublets_removed
    percent_removed = (doublets_removed / total_cells) * 100
    
    print("\n" + "="*60)
    print("SCRUBLET DOUBLET DETECTION REPORT")
    print("="*60)
    print(f"Expected doublet rate: {expected_doublet_rate:.1%}")
    print(f"Auto threshold: {scrub.threshold_:.4f}")
    print(f"Total cells: {total_cells:,}")
    print(f"Doublets removed: {doublets_removed:,} ({percent_removed:.2f}%)")
    print(f"Cells remaining: {cells_remaining:,}")
    print("="*60)
    

# 1. Load data and run initial Scrublet analysis
matrix, barcodes, features = load_data(data_dir)
scrub, doublet_scores, predicted_doublets = run_scrublet_analysis(matrix, expected_doublet_rate)
scrub.plot_histogram()
save_results(barcodes, doublet_scores, scrub.predicted_doublets_, 
              Path(output_dir) / "scrublet_results.csv")
filtered_dir = export_filtered_data(matrix, barcodes, features, scrub, 
                                   Path(output_dir) / "filtered_scrublet")
generate_final_report(matrix, scrub, expected_doublet_rate)
