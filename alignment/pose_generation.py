"""Generate scan poses for rigid and flexible axis pairs."""

from alignment.geometry import PairAligner
from alignment.polar_scan import PositionMerger, spherical_position
from molecule_utilities.logging import get_logger
from alignment import LOGGER_NAMESPACE

logger = get_logger("alignment.aligner", namespace=LOGGER_NAMESPACE)


class PoseGenerator:
    def __init__(self, reader, rigid_rows, flexible_rows, frames, scans, merge_tolerance, orientation, optimizer):
        self.reader = reader
        self.rigid_rows = rigid_rows
        self.flexible_rows = flexible_rows
        self.frames = frames
        self.scans = scans
        self.merge_tolerance = merge_tolerance
        self.orientation = orientation
        self.optimizer = optimizer

    def generate(self):
        reader, rigid_rows, flexible_rows = self.reader, self.rigid_rows, self.flexible_rows
        frames, scans, merge_tolerance = self.frames, self.scans, self.merge_tolerance
        orientation, optimizer = self.orientation, self.optimizer
        aligner = PairAligner()
        axis_pair = 0
        for r in rigid_rows:
            frame = frames.get(r.line)
            for f in flexible_rows:
                axis_pair += 1
                pose_index = 0
                prior_positions = PositionMerger(merge_tolerance) if len(scans) > 1 else None
                for scan_index, scan in enumerate(scans, 1):
                    first_plan, second_plan, third_plan = scan.plans
                    polar_merger = PositionMerger(scan.tolerance)
                    current_positions = (PositionMerger(merge_tolerance)
                                         if scan_index < len(scans) else None)
                    polar_duplicates = cross_duplicates = 0
                    for first in first_plan.distances():
                        if first is None:
                            first = (aligner.radius(r, reader.molecule(r))
                                     + aligner.radius(f, reader.molecule(f)))
                        for second in second_plan.distances():
                            for third in third_plan.distances():
                                if scan.system == "polar":
                                    point = spherical_position(first, second, third, scan.phi_unit, scan.theta_unit)
                                    distance, x, z = point.y, point.x, point.z
                                    if not polar_merger.add((x, distance, z)):
                                        polar_duplicates += 1
                                        continue
                                else:
                                    distance, x, z, point = first, second, third, None
                                position = (x, distance, z)
                                if prior_positions is not None and prior_positions.contains(position):
                                    cross_duplicates += 1
                                    continue
                                if current_positions is not None:
                                    current_positions.remember(position)
                                pose_index += 1
                                offset = None if frame is None else x * frame.x + z * frame.z
                                mol, matrix, translation, target, separation = aligner.align(
                                    r, f, reader.molecule(r), reader.molecule(f), distance, offset, orientation
                                )
                                result = None
                                if optimizer.enabled:
                                    result = optimizer.optimize(r, f, reader.molecule(r), reader.molecule(f),
                                                                mol, matrix, target)
                                    matrix = result.matrix @ matrix
                                    translation = target - matrix @ f.center
                                yield (axis_pair, pose_index, scan_index, r, f, mol, matrix, translation, target,
                                       separation, x, z, frame, orientation, point,
                                       scan.phi_unit if point is not None else "",
                                       scan.theta_unit if point is not None else "",
                                       result, optimizer.score, optimizer.objective)
                    if prior_positions is not None and current_positions is not None:
                        prior_positions.extend(current_positions)
                    if polar_duplicates:
                        logger.warning("axis pair %d: merged %d polar positions within %g Å",
                                       axis_pair, polar_duplicates, scan.tolerance)
                    if cross_duplicates:
                        logger.warning("axis pair %d scan %d: merged %d positions from earlier scans within %g Å",
                                       axis_pair, scan_index, cross_duplicates, merge_tolerance)
