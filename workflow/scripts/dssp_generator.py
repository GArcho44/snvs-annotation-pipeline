import os
import subprocess

# Function to run DSSP on a PDB file
def run_dssp(pdb_file, output_dir):
    """Run DSSP on a given PDB file."""
    # Get the base name of the PDB file (e.g., A3DHB8 from A3DHB8.pdb)
    base_name = os.path.basename(pdb_file).replace('.pdb', '')

    # Construct the output DSSP file path
    dssp_file = os.path.join(output_dir, f"{base_name}.dssp")

    # Run the mkdssp command
    try:
        subprocess.run(["mkdssp", pdb_file, dssp_file], check=True)
        print(f"Successfully generated: {dssp_file}")
    except subprocess.CalledProcessError as e:
        print(f"Failed to generate DSSP for: {base_name}.pdb")
        print(e)

# Main function to process all PDB files in the input directory
def process_all_pdb_files(input_dir, output_dir):
    """Process all PDB files in the input directory."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get all PDB files in the input directory
    pdb_files = [f for f in os.listdir(input_dir) if f.endswith('.pdb')]

    # Run DSSP on all PDB files
    for pdb_file in pdb_files:
        pdb_path = os.path.join(input_dir, pdb_file)
        run_dssp(pdb_path, output_dir)

# Main function
if __name__ == "__main__":
    # Paths
    pdb_input_dir = snakemake.input["structures_dir"]
    dssp_output_dir = snakemake.output["dssp_dir"]

    # Process all PDB files
    process_all_pdb_files(pdb_input_dir, dssp_output_dir)

    # Process complete message
    print(f"DSSP data saved to {dssp_output_dir}")
