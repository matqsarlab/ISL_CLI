import re
from typing import Optional, Any, Dict
import json
import pandas as pd
from lib.helper_functions import (
    smiles_to_mop,
    check_job_ended_norm,
    smiles_to_mol2,
    out_to_mol,
    convert_to_canonical_smiles,
    calc_mopac_desc,
)
from lib.run_software_calc import alva_desc, run_mop, ALVA_DESC_N
from pathlib import Path
from functools import cached_property
from contextlib import chdir
import logging


LOG_N = "log.log"
BASE_DIR = Path.cwd()


class ComputeDescriptors:
    def __init__(
        self,
        input_path: str,
        output_dir: str,
        out_file_name: str,
    ) -> None:
        self.input_path = Path(input_path)
        self.output_dir = Path(output_dir)
        self.output_file_name = self.output_dir / out_file_name
        self.molec_smi = "molec.smi"
        self.molec_canon_smi = "molec_canon.smi"
        self.all_data = []
        self.y_unit = str(self.smi_df.columns[0])
        self.y_file_name = self.safe_path_name(self.y_unit)

    def __call__(self) -> Any:
        print("--" * 25)
        print("\tHello from ISL!")
        print("--" * 25)
        print(f"Calculating with input: {self.input_path}")
        print(f"Output directory: {self.output_dir}")
        print(f"Output file: {self.output_file_name}")

        # Run _run_optimization only if root_dir does not exist
        if not self.output_dir.exists():
            self._run_optimization()
            print("OPTIMIZATION COMPLETED")
        if self.output_dir.is_dir() and any(
            self.output_dir.iterdir()
        ):  # Checks if directory exists and is not empty
            print("CALC_DESC")
            self._calc_descr()

    @cached_property
    def smi_df(self) -> pd.DataFrame:
        return pd.read_csv(self.input_path, index_col=0)

    def _run_optimization(self) -> None:
        df = self.smi_df
        PREDICTIONS_DIR = BASE_DIR / self.output_dir
        for i, smi in enumerate(df.index):
            canon_smi = convert_to_canonical_smiles(smi)
            print(f"\ninput SMILES: {smi}\ncanon SMILES: {canon_smi}")

            new_structure_dir = PREDICTIONS_DIR / f"{i}"
            new_structure_dir.mkdir(parents=True, exist_ok=True)

            with chdir(new_structure_dir):
                mop_file = smiles_to_mop(smi)
                out_file = run_mop(mop_file)
                end_norm_file = check_job_ended_norm(out_file)
                yval2str = str(df.iloc[i, 0])

                if end_norm_file == "SETPI":
                    file_for_alva = smiles_to_mol2(smi)
                elif end_norm_file:
                    file_for_alva = out_to_mol(end_norm_file)
                else:
                    continue  # Pomiń ten SMILES, jeśli MOPAC zawiedzie

                alva_desc(file_for_alva)

                with open(self.molec_smi, "w") as smi_name:
                    smi_name.write(smi)
                with open(self.molec_canon_smi, "w") as smi_name:
                    smi_name.write(canon_smi)
                with open(self.y_file_name, "w") as yval:
                    yval.write(yval2str)

    def _calc_descr(self) -> None:
        self.error_logger = self._setup_logger(self.output_dir / LOG_N)
        ROOT_DIR_P = Path(self.output_dir)
        PARENT_DIR_P = ROOT_DIR_P.parent
        for alva_path in ROOT_DIR_P.rglob(ALVA_DESC_N):
            molecule_record = self._process_molecule_data(alva_path)
            if molecule_record is not None:
                self.all_data.append(molecule_record)
        if self.all_data:
            self.dataframe = pd.DataFrame(self.all_data)
        else:
            self.dataframe = pd.DataFrame()
        self.dataframe.to_csv(PARENT_DIR_P / self.output_file_name, index=False)

    def _process_molecule_data(self, alva_path: Path) -> Optional[Dict[str, Any]]:
        dir_path = alva_path.parent
        alvMopac_path = dir_path / "alvaDesc_and_mopacDesc.json"
        out2_path = dir_path / "out.out"
        smi_path = dir_path / self.molec_smi
        smi_can_path = dir_path / self.molec_canon_smi
        y_path = dir_path / self.y_file_name
        try:
            molec_smi = self._read_by_logger(smi_path, "Plik SMILES")
            molec_canon_smi = self._read_by_logger(
                smi_can_path, "Plik kanonicznego SMILES"
            )
            y_val = self._read_by_logger(y_path, "Plik wartości Y")

            if molec_smi is None or molec_canon_smi is None or y_val is None:
                return None

            if not out2_path.exists():
                raise FileNotFoundError(f"Plik out.out nie znaleziony: {out2_path}")

            calc_mopac_desc(out2_path, alva_path)

            with alvMopac_path.open("r", encoding="utf-8") as file:
                descriptor_data = json.load(file)

            extracted_values = {
                desc: descriptor_data.get(desc, {}).get("value", None)
                for desc in descriptor_data
            }

            molecule_record = {
                "SMILES": molec_smi,
                "CANON_SMILES": molec_canon_smi,
                self.y_unit: y_val,
                **extracted_values,
            }
            return molecule_record

        except FileNotFoundError as fnfe:
            self.error_logger.error(f"Brak pliku w {dir_path}: {fnfe}")
            print(f"[BŁĄD] Brak wymaganego pliku w {dir_path}: {fnfe}. Pomiń.")
            return None

        except json.JSONDecodeError as jde:
            self.error_logger.exception(
                f"Błąd parsowania JSON w pliku {alvMopac_path}: {jde}"
            )
            print(f"[BŁĄD] Błąd parsowania JSON w {alvMopac_path}: {jde}. Pomiń.")
            return None

        except Exception as e:
            self.error_logger.exception(
                f"Ogólny błąd przetwarzania folderu: {dir_path}\nALVA plik: {alva_path}"
            )
            print(f"[BŁĄD] Przetwarzanie {dir_path} nie powiodło się: {e}. Pomiń.")
            return None

    def _read_by_logger(self, file_path: Path, file_description: str) -> Optional[str]:
        if not file_path.exists():
            error_msg = f"{file_description} nie znaleziony: {file_path}"
            self.error_logger.error(error_msg)
            print(f"[BŁĄD] {error_msg}. Pomiń przetwarzanie tego rekordu.")
            return None  # Zwróć None, aby zasygnalizować błąd
        return file_path.read_text(encoding="utf-8").strip()

    @staticmethod
    def safe_path_name(name: str) -> str:
        return re.sub(r'[<>:"/\\|?*\x00-\x1F]', "", name).strip()

    @staticmethod
    def _setup_logger(log_file_path: Path) -> logging.Logger:
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.ERROR)
        if not logger.handlers:
            file_handler = logging.FileHandler(log_file_path, encoding="utf-8")
            formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        return logger
