# Blind-Docking-for-Target-Fishing-based-on-Autodock
An autodock-based blind docking Target Fishing script. Used when protein loci are not known. Can be run on multiple cores.

## 1. Introduction
This is an AutoDock-based blind docking target fishing script used for docking when the protein binding sites are unknown. It can be run on multiple cores.

## 2. Principle
The implementation principle of this script is to automatically segment the receptor into multiple parts, perform docking separately, and output the score for each pose of each part.

## 3. Features
- Automatically processes multiple protein files, generates grid files, and performs docking on different regions of each protein.
- Utilizes multi-threading to accelerate the processing.
- Outputs docking results and scores to files.

## 4. Usage
1. Preprocess the ligand and save it in PDBQT format.
2. Preprocess the receptor files by breaking them down into individual chains, saving them in PDBQT format, and storing them in a folder.
3. Set the following variables in the script:
   - `ligand_file`: Path to the PDBQT file of the small molecule ligand
   - `protein_dir`: Directory path containing the protein PDBQT files
   - `output_dir`: Directory path for output results
   - `overlap`: Box overlap size (Å) (recommended to be adjusted based on the size of the small molecule ligand)
   - `cpu`: Number of cores to be used per process (recommended to be 4, adjustable based on computer configuration)
   - `exhaustiveness`: Default is 8; it is recommended to read the explanation in the official documentation and adjust as needed
   - `processes`: Total number of cores to be used

4. Run the script:
   ```bash
   python dock.py
   ```

## 5. License
MIT License
