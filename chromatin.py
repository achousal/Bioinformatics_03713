import subprocess

def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        exit(1)

def main():
    # Get file inputs from the user
    hal_file = input("Enter the path to the HAL file: ")
    genome1 = input("Enter the source genome (e.g., Human): ")
    peaks1 = input(f"Enter the path to the {genome1} peaks BED file: ")
    genome2 = input("Enter the target genome (e.g., Mouse): ")
    peaks2 = input(f"Enter the path to the {genome2} peaks BED file: ")

    output_liftover1 = f"{genome1}_to_{genome2}_liftover.bed"
    output_liftover2 = f"{genome2}_to_{genome1}_liftover.bed"
    output_overlap = "overlapping_peaks.bed"

    print("\nStarting liftover from Genome1 to Genome2...")
    run_command(f"halLiftover {hal_file} {genome1} {peaks1} {genome2} {output_liftover1}")

    print("Starting liftover from Genome2 to Genome1...")
    run_command(f"halLiftover {hal_file} {genome2} {peaks2} {genome1} {output_liftover2}")

    print("Performing intersection to find overlapping peaks...")
    run_command(f"bedtools intersect -a {output_liftover1} -b {output_liftover2} > {output_overlap}")

    print(f"Pipeline completed! Results saved in {output_overlap}")

if __name__ == "__main__":
    main()