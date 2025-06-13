from rdkit import Chem
import subprocess
import json
from datetime import date
from typing import Generator
from pathlib import Path
from lib.constants import ALVA_MOPAC_DESC_FILE_N, LOG_FILE_N


def change_keywords(
    filename: str,
    old_string: str = "PUT KEYWORDS HERE",
    keywords: str = "PM7 PRECISE PDBOUT",
) -> None:
    """
    Finds a substring within a file and replaces it with another.
    """
    with open(filename, "r") as file:
        file_contents = file.read()

    new_contents = file_contents.replace(old_string, keywords)

    with open(filename, "w") as file:
        file.write(new_contents)


def smiles_to_mop(smiles_string: str) -> str:
    """Converts a SMILES string into a MOPAC .mop file using obabel.
    As there as special characters in SMILES meaning it can break the code
    the output filename could be just same for every calc or id of calc.
    For Now same for all calc.
    Args:
        smiles_string (str): The SMILES representation of the molecule.
    Returns:
        str: The name of the generated MOPAC file.
    """
    output_filename = "out.mop"

    # Construct the obabel command
    obabel_command = (
        f'obabel -ismi -:"{smiles_string}" -omop -O "{output_filename}" --gen3d -h'
    )
    # Execute the obabel command
    subprocess.run(obabel_command, shell=True, check=True)

    change_keywords(output_filename)
    return output_filename


def check_job_ended_norm(filename: str):
    with open(filename, "r") as file:
        lines = file.readlines()

    # Check lines in reverse order
    for line in lines:
        if "AN ERROR IN ASSIGNING PI-BONDS" in line:
            return "SETPI"
        elif "Error and normal" in line:
            return False
        elif "JOB ENDED NORMALLY" in line:
            return filename

    return False  # Not found if we reach here


def smiles_to_mol2(smiles_string: str) -> str:
    base_filename = smiles_string.replace(" ", "_")
    output_filename = base_filename + ".mol2"
    obabel_command = (
        f'obabel -ismi -:"{smiles_string}" -omol2 -O "{output_filename}" --gen3d -h'
    )
    subprocess.run(obabel_command, shell=True, check=True)
    return output_filename


def out_to_mol(out_filename: str) -> str:
    mol_filename = out_filename.replace(".out", ".mol")

    obabel_command = f'obabel -iout "{out_filename}" -omol -O "{mol_filename}"'
    subprocess.run(obabel_command, shell=True, check=True)

    return mol_filename


def generate_descriptors(res_out, res_desc_names, software_version) -> Generator:
    if not res_out or not res_out[0]:
        raise ValueError("res_out is empty or not properly formatted")

    values = res_out[0]
    today = date.today().strftime("%d.%m.%Y")

    for name, value in zip(res_desc_names, values):
        yield name, {
            "value": value,
            "unit": "None",
            "metadata": {
                "software_name": "AlvaDesc",
                "software_version": software_version,
                "date": today,
            },
        }


def write_descriptors_to_json(path, generator):
    with open(path, "w") as f:
        f.write("{")
        first = True
        for name, desc in generator:
            if not first:
                f.write(",")
            else:
                first = False
            json.dump(name, f)
            f.write(":")
            json.dump(desc, f)
        f.write("}")


def convert_to_canonical_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)  # Convert SMILES to RDKit molecule object
    if mol:
        return Chem.MolToSmiles(mol, canonical=True)  # Convert back to canonical SMILES
    else:
        return smiles


def extract_float(file_lines, keyword, index):
    try:
        line = next(l for l in file_lines if keyword in l)
        return float(line.split()[index])
    except (StopIteration, IndexError, ValueError):
        return float("nan")


def calc_mopac_desc(filename: Path, json_path: str | Path):
    file_lines = [line.strip() for line in filename.read_text().splitlines()]

    hof = extract_float(file_lines, "FINAL HEAT OF FORMATION", 5)
    area = extract_float(file_lines, "COSMO AREA", 3)
    volume = extract_float(file_lines, "COSMO VOLUME", 3)
    ionisation_pot = extract_float(file_lines, "IONIZATION POTENTIAL", 3)
    homo = extract_float(file_lines, "HOMO LUMO ENERGIES (EV)", 5)
    lumo = extract_float(file_lines, "HOMO LUMO ENERGIES (EV)", 6)
    mol_weight = extract_float(file_lines, "MOLECULAR WEIGHT", 3)

    dip_x = extract_float(file_lines, "POINT-CHG.", 1)
    dip_y = extract_float(file_lines, "POINT-CHG.", 2)
    dip_z = extract_float(file_lines, "POINT-CHG.", 3)
    dip_total = extract_float(file_lines, "POINT-CHG.", 4)

    eig_idx = stop_idx = None
    for i, line in enumerate(file_lines):
        if "EIGENVALUES" in line:
            eig_idx = i
        elif "NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS" in line:
            stop_idx = i - 2
            break

    mapped = []
    if eig_idx is not None and stop_idx is not None:
        for line in file_lines[eig_idx + 1 : stop_idx]:
            for val in line.split():
                try:
                    mapped.append(float(val))
                except ValueError:
                    continue

    qmin = min(mapped) if mapped else float("nan")
    qmax = max(mapped) if mapped else float("nan")

    mopac_desc_dict = {
        "mopac_hof": hof,
        "mopac_area": area,
        "mopac_volume": volume,
        "mopac_ionisation_pot": ionisation_pot,
        "mopac_homo": homo,
        "mopac_lumo": lumo,
        "mopac_mol_weight": mol_weight,
        "mopac_dip_x": dip_x,
        "mopac_dip_y": dip_y,
        "mopac_dip_z": dip_z,
        "mopac_dip_total": dip_total,
        "mopac_qmin": qmin,
        "mopac_qmax": qmax,
    }

    detailed_dict = transform_mopac_to_detailed_dict(
        mopac_desc_dict, software_name="MOPAC", software_version="v22.1.1"
    )

    json_path = Path(json_path)
    try:
        existing_data = json.loads(json_path.read_text())
        existing_data.update(detailed_dict)
        new_path = json_path.parent / ALVA_MOPAC_DESC_FILE_N
        new_path.write_text(json.dumps(existing_data, indent=4))
        return detailed_dict
    except Exception as e:
        err_file = json_path.parents[2] / LOG_FILE_N
        with open(err_file, "a", encoding="utf-8") as f:
            f.write(f"{e}\n")
        return None


def transform_mopac_to_detailed_dict(desc_dict, software_name, software_version):
    """
    Transform a dictionary of MOPAC descriptors into a detailed structure.

    Args:
        desc_dict (dict): Dictionary containing MOPAC descriptors.
        software_name (str): Name of the software.
        software_version (str): Version of the software.

    Returns:
        dict: Transformed dictionary with metadata and value structure.
    """
    detailed_dict = {}
    current_date = date.today().strftime("%d.%m.%Y")  # Format the current date

    for key, value in desc_dict.items():
        detailed_dict[key] = {
            "value": value,
            "unit": "None",  # You can replace "None" with actual units if available
            "metadata": {
                "software_name": software_name,
                "software_version": software_version,
                "date": current_date,
            },
        }

    return detailed_dict
