import sys
import subprocess
from alvaDescCLIWrapper.alvadesccliwrapper.alvadesc import AlvaDesc
from lib.helper_functions import (
    generate_descriptors,
    write_descriptors_to_json,
)
from lib.constants import PROJECT_ROOT_P, ALVA_DESC_FILE_N


def alva_desc(filename: str):
    """ """
    match sys.platform:
        case "darwin":  # macOS
            alva_path = PROJECT_ROOT_P / "AlvaDesc_MacOS" / "alvaDescCLI"
        case "linux" | "linux2":
            alva_path = PROJECT_ROOT_P / "AlvaDesc_Linux" / "bin" / "alvaDescCLI"
        case "win32":
            alva_path = PROJECT_ROOT_P / "AlvaDesc_Windows" / "alvaDescCLI.exe"
        case _:
            raise OSError(f"Unsupported operating system: {sys.platform}")
    aDesc = AlvaDesc(str(alva_path))
    aDesc.set_input_file(filename, "MDL")
    if not aDesc.calculate_descriptors("ALL"):
        print("Error: " + aDesc.get_error())

    res_out = aDesc.get_output()
    res_desc_names = aDesc.get_output_descriptors()
    gen = generate_descriptors(res_out, res_desc_names, "alvaDesc v2.0.16")
    write_descriptors_to_json(ALVA_DESC_FILE_N, gen)
    return None


def run_mop(filename: str) -> str:
    match sys.platform:
        case "darwin":  # macOS
            mopac_path = PROJECT_ROOT_P / "mopac_mac" / "bin" / "mopac"
        case "linux" | "linux2":
            mopac_path = PROJECT_ROOT_P / "mopac_linux" / "bin" / "mopac"
        case "win32":
            mopac_path = PROJECT_ROOT_P / "mopac_windows" / "bin" / "mopac.exe"
        case _:
            raise OSError(f"Unsupported operating system: {sys.platform}")
    mopac_command = f'{mopac_path} "{filename}"'
    subprocess.run(mopac_command, shell=True, check=True)
    out_filename = filename.replace(".mop", ".out")
    return out_filename
