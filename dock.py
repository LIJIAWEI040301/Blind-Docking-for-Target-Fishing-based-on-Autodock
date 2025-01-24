import os
import subprocess
from multiprocessing import Pool
from tqdm import tqdm  
import time

# file path
ligand_file = ""  #pdbqt
protein_dir = ""  #single protein chains in pdbqt
output_dir = ""
box_size = [30.0, 30.0, 30.0]
overlap = 10.0  #box overlap size (recommended to be adjusted to the size of the small molecule ligand)
score_file = "docking_scores.txt"


os.makedirs(output_dir, exist_ok=True)


def process_region(args):
    protein_file, region_id, center_x, center_y, center_z = args
    protein_path = os.path.join(protein_dir, protein_file)
    output_file = os.path.join(output_dir, f"{protein_file}_region_{region_id}_out.pdbqt")
    log_file = os.path.join(output_dir, f"{protein_file}_region_{region_id}_log.txt")

    vina_command = [
        "vina",
        "--receptor", protein_path,
        "--ligand", ligand_file,
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(box_size[0]),
        "--size_y", str(box_size[1]),
        "--size_z", str(box_size[2]),
        "--out", output_file,
        "--log", log_file,
        "--cpu", "4",  #Number of cores to be used per process
        "--exhaustiveness", "8"  #Default is 8, it is recommended to read the explanation in the official documentation
    ]

    try:
        subprocess.run(vina_command, capture_output=True, text=True, check=True)
        scores = []
        with open(log_file, "r") as log:
            for line in log:
                if line.strip().startswith("-----+"):
                    for pose_number in range(1, 6):  # Choose the top five scoring poses
                        next_line = next(log, None)
                        if next_line:
                            try:
                                score = float(next_line.split()[1])
                                scores.append((pose_number, score))
                            except (IndexError, ValueError):
                                break
        return (protein_file, region_id, scores)
    except Exception as e:
        print(f"Error in docking {protein_file}, region {region_id}: {e}")
        return None


def save_scores_to_file(result):  #save
    if result:
        protein_file, region_id, scores = result
        with open(score_file, "a") as f:  
            for pose_number, score in scores:
                protein_name = os.path.splitext(protein_file)[0]  
                f.write(f"{protein_name} num{region_id} pose{pose_number} {score}\n")

#Protein docking site processing
def main():
    with open(score_file, "w") as f:
        f.write("name num pose score\n")  
    tasks = []

    
    for protein_file in os.listdir(protein_dir):
        if protein_file.endswith(".pdbqt"):
            protein_path = os.path.join(protein_dir, protein_file)

            
            with open(protein_path, "r") as f:
                x_coords, y_coords, z_coords = [], [], []
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        try:
                            x, y, z = map(float, line[30:54].split())
                            x_coords.append(x)
                            y_coords.append(y)
                            z_coords.append(z)
                        except ValueError:
                            continue

            if not x_coords or not y_coords or not z_coords:
                continue

            x_min, x_max = min(x_coords), max(x_coords)
            y_min, y_max = min(y_coords), max(y_coords)
            z_min, z_max = min(z_coords), max(z_coords)

            
            x_regions = int((x_max - x_min) / (box_size[0] - overlap)) + 1
            y_regions = int((y_max - y_min) / (box_size[1] - overlap)) + 1
            z_regions = int((z_max - z_min) / (box_size[2] - overlap)) + 1

            region_id = 1
            for i in range(x_regions):
                for j in range(y_regions):
                    for k in range(z_regions):
                        center_x = x_min + i * (box_size[0] - overlap) + box_size[0] / 2
                        center_y = y_min + j * (box_size[1] - overlap) + box_size[1] / 2
                        center_z = z_min + k * (box_size[2] - overlap) + box_size[2] / 2
                        tasks.append((protein_file, region_id, center_x, center_y, center_z))
                        region_id += 1

    # multiprocess
    start_time = time.time()  
    with Pool(processes=64) as pool:  # Total cores to be used
        # progress bar
        for result in tqdm(pool.imap(process_region, tasks), total=len(tasks), desc="Processing regions"):
            save_scores_to_file(result)  

    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Docking completed in {elapsed_time:.2f} seconds.")


if __name__ == "__main__":
    main()
