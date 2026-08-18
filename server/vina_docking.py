"""Blind cholesterol docking and collision-free PDB assembly for CholBindNet.

The workflow follows the supplied AutoDock Vina notebook:

1. prepare a heterogen-free receptor copy for Vina,
2. prepare the receptor with ADFRsuite,
3. generate a cholesterol conformer with RDKit and prepare it with Meeko,
4. create a whole-protein docking box with 10 Angstrom padding,
5. run Vina with exhaustiveness 16 and up to 100 modes, and
6. append the docked poses as correctly formatted CLR residues to the untouched
   uploaded PDB, preserving every original experimental ligand and residue name.

External tool paths can be configured with environment variables documented by
``DockingConfig.from_environment``.  All subprocesses run without a shell and
inside a per-request working directory.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence


CHOLESTEROL_SMILES = (
    "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4"
    "[C@@]3(CC[C@@H](C4)O)C)C"
)

PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "ASX", "GLX", "SEC", "PYL", "MSE",
}

CHAIN_IDS = "LXYZABCDEFGHIJKLMNOPQRSTUVWabcdefghijklmnopqrstuvwxyz0123456789"
VINA_RESULT_RE = re.compile(
    r"REMARK\s+VINA\s+RESULT:\s*([-+]?\d+(?:\.\d+)?)",
    re.IGNORECASE,
)


class DockingError(RuntimeError):
    """Raised when receptor preparation or docking cannot be completed."""


def _env_int(name: str, default: int) -> int:
    value = os.environ.get(name)
    if value is None:
        return default
    try:
        return int(value)
    except ValueError as exc:
        raise DockingError(f"{name} must be an integer; received {value!r}.") from exc


def _env_float(name: str, default: float) -> float:
    value = os.environ.get(name)
    if value is None:
        return default
    try:
        return float(value)
    except ValueError as exc:
        raise DockingError(f"{name} must be numeric; received {value!r}.") from exc


@dataclass(frozen=True)
class DockingConfig:
    """Runtime configuration for the notebook-derived Vina workflow.

    Environment variables:

    - ``ADFR_PREPARE_RECEPTOR``: ADFRsuite ``prepare_receptor`` executable.
    - ``MEEKO_PREPARE_LIGAND``: Meeko ``mk_prepare_ligand.py`` script.
    - ``VINA_BIN``: AutoDock Vina executable.
    - ``CHOLNET_VINA_EXHAUSTIVENESS``: defaults to 16.
    - ``CHOLNET_VINA_NUM_MODES``: defaults to 100.
    - ``CHOLNET_VINA_PADDING``: whole-protein box padding; defaults to 10 A.
    - ``CHOLNET_VINA_TIMEOUT``: per-command timeout; defaults to 1800 seconds.
    - ``CHOLNET_VINA_CHAIN``: preferred ligand chain; defaults to L.
    """

    prepare_receptor: Optional[str] = None
    prepare_ligand: Optional[str] = None
    vina: Optional[str] = None
    exhaustiveness: int = 16
    num_modes: int = 100
    padding: float = 10.0
    timeout_seconds: int = 1800
    preferred_chain: str = "L"

    @classmethod
    def from_environment(cls, base_dir: Optional[str] = None) -> "DockingConfig":
        base = Path(base_dir or Path(__file__).resolve().parent)
        default_prepare_ligand = base / "mk_prepare_ligand.py"
        preferred_chain = os.environ.get("CHOLNET_VINA_CHAIN", "L").strip() or "L"
        if len(preferred_chain) != 1 or not preferred_chain.isalnum():
            raise DockingError("CHOLNET_VINA_CHAIN must be one alphanumeric character.")
        return cls(
            prepare_receptor=os.environ.get("ADFR_PREPARE_RECEPTOR"),
            prepare_ligand=os.environ.get(
                "MEEKO_PREPARE_LIGAND",
                str(default_prepare_ligand) if default_prepare_ligand.exists() else None,
            ),
            vina=os.environ.get("VINA_BIN"),
            exhaustiveness=_env_int("CHOLNET_VINA_EXHAUSTIVENESS", 16),
            num_modes=_env_int("CHOLNET_VINA_NUM_MODES", 100),
            padding=_env_float("CHOLNET_VINA_PADDING", 10.0),
            timeout_seconds=_env_int("CHOLNET_VINA_TIMEOUT", 1800),
            preferred_chain=preferred_chain,
        )


@dataclass(frozen=True)
class PoseAtom:
    name: str
    x: float
    y: float
    z: float
    element: str


@dataclass(frozen=True)
class VinaPose:
    rank: int
    atoms: tuple[PoseAtom, ...]
    affinity_kcal_mol: Optional[float] = None


@dataclass(frozen=True)
class AssignedPose:
    rank: int
    chain_id: str
    residue_number: int
    atom_count: int
    affinity_kcal_mol: Optional[float] = None

    def as_dict(self) -> dict:
        return {
            "rank": self.rank,
            "chain_id": self.chain_id,
            "residue_number": self.residue_number,
            "atom_count": self.atom_count,
            "affinity_kcal_mol": self.affinity_kcal_mol,
        }


@dataclass(frozen=True)
class DockingResult:
    complex_pdb: str
    all_poses_pdb: str
    cleaned_receptor_pdb: str
    receptor_pdbqt: str
    vina_output_pdbqt: str
    assigned_poses: tuple[AssignedPose, ...]
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    config: DockingConfig
    warnings: tuple[str, ...] = ()

    def api_metadata(self) -> dict:
        return {
            "method": "AutoDock Vina blind docking",
            "pose_count": len(self.assigned_poses),
            "box": {
                "center": {
                    "x": self.center[0],
                    "y": self.center[1],
                    "z": self.center[2],
                },
                "size": {
                    "x": self.size[0],
                    "y": self.size[1],
                    "z": self.size[2],
                },
                "padding_angstrom": self.config.padding,
            },
            "parameters": {
                "exhaustiveness": self.config.exhaustiveness,
                "requested_num_modes": self.config.num_modes,
            },
            "poses": [pose.as_dict() for pose in self.assigned_poses],
            "warnings": list(self.warnings),
        }


def _resolve_executable(configured: Optional[str], name: str, fallbacks: Sequence[str]) -> str:
    candidates = [configured, shutil.which(name), *fallbacks]
    for candidate in candidates:
        if not candidate:
            continue
        expanded = os.path.abspath(os.path.expanduser(candidate)) if os.path.sep in candidate else candidate
        resolved = shutil.which(expanded) or (expanded if os.path.isfile(expanded) else None)
        if resolved and os.access(resolved, os.X_OK):
            return resolved
    variable = {
        "prepare_receptor": "ADFR_PREPARE_RECEPTOR",
        "vina": "VINA_BIN",
    }.get(name, name.upper())
    raise DockingError(
        f"Required executable {name!r} was not found. Set {variable} to its full path."
    )


def _resolve_script(configured: Optional[str], name: str, fallbacks: Sequence[str]) -> str:
    candidates = [configured, shutil.which(name), *fallbacks]
    for candidate in candidates:
        if not candidate:
            continue
        expanded = os.path.abspath(os.path.expanduser(candidate)) if os.path.sep in candidate else candidate
        resolved = shutil.which(expanded) or (expanded if os.path.isfile(expanded) else None)
        if resolved and os.path.isfile(resolved):
            return resolved
    raise DockingError(
        "Meeko's mk_prepare_ligand.py was not found. Set MEEKO_PREPARE_LIGAND "
        "to its full path."
    )


def _run_command(
    command: Sequence[str],
    *,
    cwd: str,
    timeout_seconds: int,
    log_path: Optional[str] = None,
) -> subprocess.CompletedProcess:
    try:
        result = subprocess.run(
            list(command),
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise DockingError(
            f"Command {os.path.basename(command[0])!r} exceeded the "
            f"{timeout_seconds}-second timeout."
        ) from exc
    except OSError as exc:
        raise DockingError(f"Could not start {command[0]!r}: {exc}") from exc

    if log_path:
        with open(log_path, "w", encoding="utf-8") as handle:
            handle.write(result.stdout or "")
            if result.stderr:
                handle.write("\n--- STDERR ---\n")
                handle.write(result.stderr)

    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "No diagnostic output").strip()
        if len(detail) > 2000:
            detail = detail[-2000:]
        raise DockingError(
            f"{os.path.basename(command[0])} failed with exit code "
            f"{result.returncode}: {detail}"
        )
    return result


def _valid_xyz(line: str) -> bool:
    try:
        float(line[30:38])
        float(line[38:46])
        float(line[46:54])
        return True
    except (ValueError, IndexError):
        return False


def _fallback_clean_receptor(source: str, destination: str) -> None:
    kept = []
    with open(source, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line.startswith("ATOM  ") or len(line) < 54 or not _valid_xyz(line):
                continue
            residue_name = line[17:20].strip().upper()
            alternate = line[16:17]
            if residue_name not in PROTEIN_RESIDUES or alternate not in {" ", "A"}:
                continue
            if alternate == "A":
                line = line[:16] + " " + line[17:]
            kept.append(line)
    if not kept:
        raise DockingError(
            "The uploaded file contains no readable standard-protein ATOM records."
        )
    with open(destination, "w", encoding="utf-8") as handle:
        handle.write("\n".join(kept))
        handle.write("\nTER\nEND\n")


def clean_receptor(source: str, destination: str) -> tuple[str, ...]:
    """Remove heterogens, with a conservative ATOM-only fallback."""
    warnings = []
    try:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer

        fixer = PDBFixer(filename=source)
        fixer.removeHeterogens(keepWater=False)
        with open(destination, "w", encoding="utf-8") as handle:
            PDBFile.writeFile(fixer.topology, fixer.positions, handle)
        if not _coordinates_from_pdb(destination):
            raise DockingError("PDBFixer produced a receptor with no atoms.")
    except (ImportError, ModuleNotFoundError) as exc:
        warnings.append(
            f"PDBFixer/OpenMM was unavailable ({exc}); used conservative ATOM-only cleaning."
        )
        _fallback_clean_receptor(source, destination)
    except Exception as exc:
        warnings.append(
            f"PDBFixer could not read the upload ({exc}); used conservative ATOM-only cleaning."
        )
        _fallback_clean_receptor(source, destination)
    return tuple(warnings)


def _coordinates_from_pdb(path: str) -> list[tuple[float, float, float]]:
    coordinates = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")) and _valid_xyz(line):
                coordinates.append(
                    (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                )
    return coordinates


def blind_docking_box(
    receptor_pdb: str, padding: float
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    coordinates = _coordinates_from_pdb(receptor_pdb)
    if not coordinates:
        raise DockingError("The cleaned receptor contains no readable coordinates.")
    mins = tuple(min(values) for values in zip(*coordinates))
    maxs = tuple(max(values) for values in zip(*coordinates))
    center = tuple((low + high) / 2.0 for low, high in zip(mins, maxs))
    size = tuple((high - low) + 2.0 * padding for low, high in zip(mins, maxs))
    if any(value <= 0 for value in size):
        raise DockingError("The computed blind-docking box has an invalid size.")
    return center, size


def _write_cholesterol_sdf(path: str) -> None:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except (ImportError, ModuleNotFoundError) as exc:
        raise DockingError(
            "RDKit is required to generate the cholesterol ligand."
        ) from exc

    molecule = Chem.MolFromSmiles(CHOLESTEROL_SMILES)
    if molecule is None:
        raise DockingError("RDKit could not parse the cholesterol SMILES string.")
    molecule = Chem.AddHs(molecule)
    status = AllChem.EmbedMolecule(molecule, randomSeed=42)
    if status != 0:
        raise DockingError("RDKit could not generate a 3D cholesterol conformer.")
    writer = Chem.SDWriter(path)
    try:
        writer.write(molecule)
    finally:
        writer.close()
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise DockingError("RDKit did not create the cholesterol SDF file.")


def _element_from_pdbqt(line: str, atom_name: str) -> str:
    fields = line.split()
    autodock_type = fields[-1].upper() if fields else ""
    type_map = {
        "A": "C", "C": "C", "OA": "O", "O": "O", "HD": "H",
        "H": "H", "N": "N", "NA": "N", "S": "S", "SA": "S",
        "P": "P", "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
    }
    if autodock_type in type_map:
        return type_map[autodock_type]
    cleaned = re.sub(r"[^A-Za-z]", "", atom_name)
    if not cleaned:
        return "C"
    if cleaned[:2].upper() in {"CL", "BR"}:
        return cleaned[:2].title()
    return cleaned[0].upper()


def _pose_atom_from_pdbqt(line: str) -> PoseAtom:
    name = line[12:16].strip() if len(line) >= 16 else ""
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except (ValueError, IndexError):
        fields = line.split()
        if len(fields) < 9:
            raise DockingError(f"Could not parse a Vina atom record: {line!r}")
        name = name or fields[2]
        try:
            x, y, z = map(float, fields[5:8])
        except ValueError as exc:
            raise DockingError(f"Could not parse Vina coordinates: {line!r}") from exc
    name = name or "C"
    return PoseAtom(name=name, x=x, y=y, z=z, element=_element_from_pdbqt(line, name))


def read_vina_poses(path: str) -> tuple[VinaPose, ...]:
    """Read heavy-atom poses and Vina affinities directly from multi-model PDBQT."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise DockingError("Vina did not create a non-empty pose file.")

    poses = []
    atoms: list[PoseAtom] = []
    affinity: Optional[float] = None

    def finish_pose() -> None:
        nonlocal atoms, affinity
        heavy_atoms = tuple(atom for atom in atoms if atom.element.upper() != "H")
        if heavy_atoms:
            poses.append(
                VinaPose(
                    rank=len(poses) + 1,
                    atoms=heavy_atoms,
                    affinity_kcal_mol=affinity,
                )
            )
        atoms = []
        affinity = None

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith("MODEL"):
                if atoms:
                    finish_pose()
                continue
            match = VINA_RESULT_RE.search(line)
            if match:
                affinity = float(match.group(1))
                continue
            if line.startswith(("ATOM  ", "HETATM")):
                atoms.append(_pose_atom_from_pdbqt(line))
                continue
            if line.startswith("ENDMDL"):
                finish_pose()
    if atoms:
        finish_pose()
    if not poses:
        raise DockingError("No atom-containing poses could be read from Vina output.")
    return tuple(poses)


def read_pose_pdb(path: str) -> tuple[VinaPose, ...]:
    """Read pose blocks separated by TER, including the notebook's malformed PDB."""
    poses = []
    atoms: list[PoseAtom] = []

    def finish_pose() -> None:
        nonlocal atoms
        if atoms:
            poses.append(VinaPose(rank=len(poses) + 1, atoms=tuple(atoms)))
        atoms = []

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith(("ATOM  ", "HETATM")):
                if not _valid_xyz(line):
                    raise DockingError(f"Could not parse pose coordinates: {line!r}")
                name = line[12:16].strip() or "C"
                element = line[76:78].strip() if len(line) >= 78 else ""
                if not element:
                    element = _element_from_pdbqt(line, name)
                if element.upper() != "H":
                    atoms.append(
                        PoseAtom(
                            name=name,
                            x=float(line[30:38]),
                            y=float(line[38:46]),
                            z=float(line[46:54]),
                            element=element,
                        )
                    )
            elif line.startswith(("TER", "ENDMDL")):
                finish_pose()
    finish_pose()
    if not poses:
        raise DockingError("No poses were found in the supplied pose PDB.")
    return tuple(poses)


def _receptor_identifiers(path: Optional[str]) -> tuple[int, set[str], set[tuple[str, int]]]:
    maximum_serial = 0
    chains: set[str] = set()
    residue_ids: set[tuple[str, int]] = set()
    if not path:
        return maximum_serial, chains, residue_ids
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                maximum_serial = max(maximum_serial, int(line[6:11]))
            except ValueError:
                pass
            chain = line[21:22].strip()
            if chain:
                chains.add(chain)
            try:
                residue_ids.add((chain, int(line[22:26])))
            except ValueError:
                pass
    return maximum_serial, chains, residue_ids


def _choose_chain(preferred: str, receptor_chains: set[str]) -> str:
    candidates = preferred + "".join(chain for chain in CHAIN_IDS if chain != preferred)
    for chain in candidates:
        if chain not in receptor_chains:
            return chain
    return preferred


def _free_residue_numbers(
    chain_id: str, used: set[tuple[str, int]], count: int
) -> list[int]:
    free = [number for number in range(1, 10000) if (chain_id, number) not in used]
    if len(free) < count:
        raise DockingError(
            f"Chain {chain_id!r} does not have {count} free four-digit residue numbers."
        )
    return free[:count]


def _format_atom_name(name: str) -> str:
    cleaned = re.sub(r"\s+", "", name)[:4] or "C"
    return f"{cleaned:<4}"


def _format_pdb_atom(
    serial: int,
    atom: PoseAtom,
    chain_id: str,
    residue_number: int,
) -> str:
    if not 1 <= serial <= 99999:
        raise DockingError(
            "Appending docked poses would exceed the PDB atom-serial limit of 99999."
        )
    return (
        f"HETATM{serial:5d} {_format_atom_name(atom.name)} CLR {chain_id}"
        f"{residue_number:4d}    {atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
        f"{1.00:6.2f}{0.00:6.2f}          {atom.element:>2}  "
    )


def write_assigned_poses(
    poses: Sequence[VinaPose],
    output_path: str,
    *,
    receptor_pdb: Optional[str] = None,
    preferred_chain: str = "L",
) -> tuple[AssignedPose, ...]:
    """Write CLR poses with identifiers that cannot clash with the receptor."""
    maximum_serial, receptor_chains, used_residues = _receptor_identifiers(receptor_pdb)
    chain_id = _choose_chain(preferred_chain, receptor_chains)
    residue_numbers = _free_residue_numbers(chain_id, used_residues, len(poses))
    serial = maximum_serial + 1
    assigned = []

    with open(output_path, "w", encoding="utf-8") as handle:
        for pose, residue_number in zip(poses, residue_numbers):
            for atom in pose.atoms:
                handle.write(
                    _format_pdb_atom(serial, atom, chain_id, residue_number) + "\n"
                )
                serial += 1
            handle.write("TER\n")
            assigned.append(
                AssignedPose(
                    rank=pose.rank,
                    chain_id=chain_id,
                    residue_number=residue_number,
                    atom_count=len(pose.atoms),
                    affinity_kcal_mol=pose.affinity_kcal_mol,
                )
            )
        handle.write("END\n")
    return tuple(assigned)


def append_poses_to_receptor(
    source_pdb: str, poses_pdb: str, output_path: str
) -> None:
    """Append docked CLR poses while preserving the uploaded PDB verbatim.

    ``source_pdb`` must be the original upload, not the heterogen-free receptor
    prepared for Vina.  New coordinate records are inserted before the PDB
    connectivity/footer section so original HETATM residue names and CONECT
    records remain intact.
    """
    source_records = []
    has_atom_records = False
    with open(source_pdb, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line.rstrip("\r\n")
            source_records.append(record)
            if record.startswith("ATOM  "):
                has_atom_records = True
    if not has_atom_records:
        raise DockingError("The uploaded PDB has no protein ATOM records to append to.")

    pose_records = []
    with open(poses_pdb, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM", "TER")):
                pose_records.append(line.rstrip("\r\n"))
    if not any(line.startswith("HETATM") for line in pose_records):
        raise DockingError("The pose PDB contains no docked ligand atoms.")

    insertion_index = len(source_records)
    for index, line in enumerate(source_records):
        if line.startswith(("CONECT", "MASTER", "ENDMDL", "END   ")) or line == "END":
            insertion_index = index
            break

    before_footer = source_records[:insertion_index]
    footer = source_records[insertion_index:]
    while before_footer and not before_footer[-1]:
        before_footer.pop()

    output_records = [
        *before_footer,
        "REMARK 900 DOCKED CLR POSES APPENDED; ORIGINAL HETATM NAMES PRESERVED",
        "TER",
        *pose_records,
        *footer,
    ]
    if not any(line == "END" or line.startswith("END   ") for line in footer):
        output_records.append("END")

    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(output_records))
        handle.write("\n")


def normalize_pose_pdb(
    input_pose_pdb: str,
    output_pose_pdb: str,
    *,
    receptor_pdb: Optional[str] = None,
    preferred_chain: str = "L",
) -> tuple[AssignedPose, ...]:
    """Correct an existing all_poses.pdb, including the supplied notebook output."""
    poses = read_pose_pdb(input_pose_pdb)
    return write_assigned_poses(
        poses,
        output_pose_pdb,
        receptor_pdb=receptor_pdb,
        preferred_chain=preferred_chain,
    )


def dock_cholesterol(
    receptor_pdb: str,
    workdir: str,
    config: Optional[DockingConfig] = None,
) -> DockingResult:
    """Run blind cholesterol docking and return a prediction-ready complex PDB."""
    receptor_pdb = os.path.abspath(receptor_pdb)
    workdir = os.path.abspath(workdir)
    if not os.path.isfile(receptor_pdb):
        raise DockingError(f"Uploaded receptor was not found: {receptor_pdb}")
    if os.path.getsize(receptor_pdb) == 0:
        raise DockingError("The uploaded receptor file is empty.")

    config = config or DockingConfig.from_environment()
    if config.exhaustiveness < 1 or config.num_modes < 1:
        raise DockingError("Vina exhaustiveness and num_modes must both be positive.")
    if config.padding < 0 or config.timeout_seconds < 1:
        raise DockingError("Vina padding must be nonnegative and timeout must be positive.")

    os.makedirs(workdir, exist_ok=True)
    prepare_receptor = _resolve_executable(
        config.prepare_receptor,
        "prepare_receptor",
        ("/shared/apps/ADFR/1.0/bin/prepare_receptor",),
    )
    prepare_ligand = _resolve_script(
        config.prepare_ligand,
        "mk_prepare_ligand.py",
        (str(Path(__file__).resolve().parent / "mk_prepare_ligand.py"),),
    )
    vina = _resolve_executable(
        config.vina,
        "vina",
        (
            "/home/aashish/software/docking/AutoDock-Vina-1.2.5/"
            "binary_code/vina_1.2.3_linux_x86_64",
        ),
    )

    cleaned_pdb = os.path.join(workdir, "receptor_clean.pdb")
    receptor_pdbqt = os.path.join(workdir, "receptor.pdbqt")
    ligand_sdf = os.path.join(workdir, "CLR.sdf")
    ligand_pdbqt = os.path.join(workdir, "CLR.pdbqt")
    vina_output = os.path.join(workdir, "vina_out.pdbqt")
    vina_log = os.path.join(workdir, "vina_log.txt")
    all_poses_pdb = os.path.join(workdir, "all_poses.pdb")
    complex_pdb = os.path.join(workdir, "complex_all_poses.pdb")

    warnings = list(clean_receptor(receptor_pdb, cleaned_pdb))
    _run_command(
        [
            prepare_receptor,
            "-r", cleaned_pdb,
            "-o", receptor_pdbqt,
            "-A", "checkhydrogens",
            "-U", "nphs_lps_waters",
        ],
        cwd=workdir,
        timeout_seconds=config.timeout_seconds,
    )

    _write_cholesterol_sdf(ligand_sdf)
    ligand_command = [prepare_ligand, "-i", ligand_sdf, "-o", ligand_pdbqt]
    if prepare_ligand.lower().endswith(".py"):
        ligand_command.insert(0, sys.executable)
    _run_command(
        ligand_command,
        cwd=workdir,
        timeout_seconds=config.timeout_seconds,
    )

    center, size = blind_docking_box(cleaned_pdb, config.padding)
    _run_command(
        [
            vina,
            "--receptor", receptor_pdbqt,
            "--ligand", ligand_pdbqt,
            "--center_x", f"{center[0]:.6f}",
            "--center_y", f"{center[1]:.6f}",
            "--center_z", f"{center[2]:.6f}",
            "--size_x", f"{size[0]:.6f}",
            "--size_y", f"{size[1]:.6f}",
            "--size_z", f"{size[2]:.6f}",
            "--exhaustiveness", str(config.exhaustiveness),
            "--num_modes", str(config.num_modes),
            "--out", vina_output,
        ],
        cwd=workdir,
        timeout_seconds=config.timeout_seconds,
        log_path=vina_log,
    )

    poses = read_vina_poses(vina_output)
    if any(len(pose.atoms) != 28 for pose in poses):
        warnings.append(
            "At least one docked pose did not contain the expected 28 cholesterol heavy atoms."
        )
    assigned = write_assigned_poses(
        poses,
        all_poses_pdb,
        # Allocate chain/residue/serial IDs against every original ATOM and
        # HETATM record so the appended docked CLR poses cannot collide with an
        # experimental ligand that was present in the uploaded structure.
        receptor_pdb=receptor_pdb,
        preferred_chain=config.preferred_chain,
    )
    # Vina must use the cleaned receptor, but the prediction/display complex
    # must use the untouched upload so experimental ligands are not discarded.
    append_poses_to_receptor(receptor_pdb, all_poses_pdb, complex_pdb)

    return DockingResult(
        complex_pdb=complex_pdb,
        all_poses_pdb=all_poses_pdb,
        cleaned_receptor_pdb=cleaned_pdb,
        receptor_pdbqt=receptor_pdbqt,
        vina_output_pdbqt=vina_output,
        assigned_poses=assigned,
        center=center,
        size=size,
        config=config,
        warnings=tuple(warnings),
    )