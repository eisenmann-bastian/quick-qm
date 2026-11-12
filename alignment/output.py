"""Transactional output of aligned SDF records and CSV metadata."""

import csv
import json
import os
import shutil
import tempfile
from pathlib import Path

from rdkit import Chem

from molecule_utilities.logging import get_logger
from alignment import LOGGER_NAMESPACE

logger = get_logger("alignment.aligner", namespace=LOGGER_NAMESPACE)


def write_output(output, output_format, pairs, destinations, overwrite):
    fields = ("pair_id", "axis_pair_id", "pose_index", "scan_index", "sdf_file", "sdf_record_index", "rigid_csv_line", "flexible_csv_line",
              "rigid_source", "flexible_source", "distance_angstrom", "rigid_center",
              "flexible_original_center", "flexible_aligned_center", "rigid_axis",
              "flexible_original_axis", "flexible_aligned_axis", "rotation", "translation",
              "x_offset_angstrom", "z_offset_angstrom", "x_axis", "z_axis",
              "reference_method", "reference_atom_index", "reference_point", "flexible_orientation",
              "coordinate_system", "polar_radius_angstrom", "polar_phi", "polar_phi_unit", "polar_theta", "polar_theta_unit",
              "torsion_angle_degrees", "torsion_score", "torsion_score_type", "torsion_objective")
    staging = Path(tempfile.mkdtemp(prefix=".aligner-", dir=output.parent))
    preserve_staging = False
    writer = None
    try:
        with (staging / "metadata.csv").open("x", encoding="utf-8", newline="") as metadata:
            if output_format == "single":
                writer = Chem.SDWriter(str(staging / "single.sdf"))
            current_axis_pair = None
            table = csv.DictWriter(metadata, fieldnames=fields, lineterminator="\n")
            table.writeheader()
            index = 0
            for index, (axis_pair, pose_index, scan_index, rigid, flexible, molecule, matrix, translation, target,
                        distance, x, z, frame, orientation, point, phi_unit, theta_unit,
                        optimization_result, score_type, objective) in enumerate(pairs, 1):
                pair_id = f"pair-{index:06d}"
                axis_pair_id = f"axis-pair-{axis_pair:06d}"
                sdf = output / f"{axis_pair_id}.sdf" if output_format == "pair-files" else output
                molecule.SetProp("_Name", pair_id)
                molecule.SetProp("quick_qm_pair_id", pair_id)
                molecule.SetProp("quick_qm_axis_pair_id", axis_pair_id)
                molecule.SetProp("quick_qm_pose_index", str(pose_index))
                molecule.SetProp("quick_qm_scan_index", str(scan_index))
                molecule.SetProp("quick_qm_distance_angstrom", str(distance))
                molecule.SetProp("quick_qm_flexible_orientation", orientation)
                molecule.SetProp("quick_qm_coordinate_system", "polar" if point is not None else "grid")
                if optimization_result is not None:
                    molecule.SetProp("quick_qm_torsion_angle_degrees", str(optimization_result.angle))
                    molecule.SetProp("quick_qm_torsion_score", str(optimization_result.score))
                    molecule.SetProp("quick_qm_torsion_score_type", score_type)
                    molecule.SetProp("quick_qm_torsion_objective", objective)
                if point is not None:
                    molecule.SetProp("quick_qm_polar_radius_angstrom", str(point.radius))
                    molecule.SetProp("quick_qm_polar_phi", str(point.phi))
                    molecule.SetProp("quick_qm_polar_phi_unit", phi_unit)
                    molecule.SetProp("quick_qm_polar_theta", str(point.theta))
                    molecule.SetProp("quick_qm_polar_theta_unit", theta_unit)
                molecule.SetProp("quick_qm_x_offset_angstrom", str(x))
                molecule.SetProp("quick_qm_z_offset_angstrom", str(z))
                if frame is not None:
                    molecule.SetProp("quick_qm_reference_method", frame.method)
                    molecule.SetProp("quick_qm_x_axis", json.dumps(frame.x.tolist()))
                    molecule.SetProp("quick_qm_z_axis", json.dumps(frame.z.tolist()))
                molecule.SetProp("quick_qm_rigid_source", json.dumps(rigid.fields, separators=(",", ":")))
                molecule.SetProp("quick_qm_flexible_source", json.dumps(flexible.fields, separators=(",", ":")))
                if output_format == "pair-files" and current_axis_pair != axis_pair:
                    if writer is not None:
                        writer.flush()
                        writer.close()
                    writer = Chem.SDWriter(str(staging / f"{axis_pair_id}.sdf"))
                    current_axis_pair = axis_pair
                writer.write(molecule)
                encode = lambda value: json.dumps(value, separators=(",", ":"), allow_nan=False)
                table.writerow({
                    "pair_id": pair_id, "axis_pair_id": axis_pair_id, "pose_index": pose_index, "scan_index": scan_index,
                    "sdf_file": os.path.relpath(sdf),
                    "sdf_record_index": pose_index - 1 if output_format == "pair-files" else index - 1,
                    "rigid_csv_line": rigid.line, "flexible_csv_line": flexible.line,
                    "rigid_source": encode(rigid.fields), "flexible_source": encode(flexible.fields),
                    "distance_angstrom": distance, "rigid_center": encode(rigid.center.tolist()),
                    "flexible_original_center": encode(flexible.center.tolist()),
                    "flexible_aligned_center": encode(target.tolist()), "rigid_axis": encode(rigid.axis.tolist()),
                    "flexible_original_axis": encode(flexible.axis.tolist()),
                    "flexible_aligned_axis": encode((matrix @ flexible.axis).tolist()),
                    "rotation": encode(matrix.tolist()), "translation": encode(translation.tolist()),
                    "x_offset_angstrom": x, "z_offset_angstrom": z,
                    "x_axis": encode(frame.x.tolist()) if frame is not None else "",
                    "z_axis": encode(frame.z.tolist()) if frame is not None else "",
                    "reference_method": frame.method if frame is not None else "",
                    "reference_atom_index": frame.atom_index if frame is not None and frame.atom_index is not None else "",
                    "reference_point": encode(frame.point.tolist()) if frame is not None else "",
                    "flexible_orientation": orientation,
                    "coordinate_system": "polar" if point is not None else "grid",
                    "polar_radius_angstrom": point.radius if point is not None else "",
                    "polar_phi": point.phi if point is not None else "",
                    "polar_phi_unit": phi_unit,
                    "polar_theta": point.theta if point is not None else "",
                    "polar_theta_unit": theta_unit,
                    "torsion_angle_degrees": optimization_result.angle if optimization_result is not None else "",
                    "torsion_score": optimization_result.score if optimization_result is not None else "",
                    "torsion_score_type": score_type if optimization_result is not None else "",
                    "torsion_objective": objective if optimization_result is not None else "",
                })
            if writer is not None:
                writer.flush()
                writer.close()
                writer = None
        staged_files = ([staging / path.name for path in destinations[:-1]]
                        if output_format == "pair-files" else [staging / "single.sdf"])
        staged_files.append(staging / "metadata.csv")
        installed, backups = [], []
        created_directory = False
        try:
            if output_format == "pair-files":
                if not overwrite or not output.exists():
                    output.mkdir()
                    created_directory = True
                elif output.is_symlink() or not output.is_dir():
                    raise ValueError("pair-files output must be a directory, not a file or symlink")
            for file_index, (source, destination) in enumerate(zip(staged_files, destinations)):
                if destination.is_dir():
                    raise ValueError(f"output file is a directory: {destination}")
                if overwrite and os.path.lexists(destination):
                    backup = staging / f"backup-{file_index}"
                    os.replace(destination, backup)
                    backups.append((backup, destination))
                if overwrite:
                    os.replace(source, destination)
                else:
                    os.link(source, destination)
                installed.append(destination)
        except (OSError, ValueError, RuntimeError):
            for path in reversed(installed):
                try:
                    path.unlink()
                except OSError as cleanup_error:
                    logger.error("could not remove incomplete output %s: %s", path, cleanup_error)
            for backup, destination in reversed(backups):
                try:
                    os.replace(backup, destination)
                except OSError as cleanup_error:
                    preserve_staging = True
                    logger.error("could not restore %s: %s; backup retained at %s",
                                 destination, cleanup_error, backup)
            if created_directory:
                try:
                    output.rmdir()
                except OSError as cleanup_error:
                    logger.error("could not remove output directory %s: %s", output, cleanup_error)
            raise
        return index
    finally:
        if writer is not None:
            writer.close()
        if not preserve_staging:
            try:
                shutil.rmtree(staging)
            except OSError as cleanup_error:
                logger.error("could not clean up temporary output %s: %s", staging, cleanup_error)
