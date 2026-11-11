"""Coordinate functional-group selection and alignment output."""

import io

from .alignment import AlignmentCsv, AlignmentSdf
from .batch import BatchSelector
from .matching import MoleculeSelector
from .output import CsvResultWriter


class SelectionPipeline:
    """Run selection and optionally derive alignment results in memory."""

    def __init__(self, molecule_selector=None):
        self.matcher = molecule_selector if molecule_selector is not None else MoleculeSelector()
        self.batch = BatchSelector(self.matcher)
        self.selection_writer = CsvResultWriter()

    def select_csv(self, path, *, delimiter=","):
        return self.batch.select(path, delimiter=delimiter)

    def prepare(self, selections, output_format, *, sdf_file=None, axis_config=None):
        if output_format == "selection":
            return selections
        if output_format not in ("alignment", "sdf"):
            raise ValueError(f"unknown output format: {output_format}")
        stream = io.StringIO(newline="")
        self.selection_writer.write(selections, stream, sdf_file=sdf_file)
        stream.seek(0)
        alignment = AlignmentCsv(axis_config)
        axes = alignment.read_stream(stream, include_unmatched=output_format == "sdf")
        if output_format == "sdf":
            return AlignmentSdf.build(alignment.records, axes)
        return axes

    def write(self, result, destination, output_format, *, sdf_file=None):
        if output_format == "selection":
            self.selection_writer.write(result, destination, sdf_file=sdf_file)
        elif output_format == "alignment":
            AlignmentCsv.write(result, destination)
        elif output_format == "sdf":
            AlignmentSdf.write(result, destination)
        else:
            raise ValueError(f"unknown output format: {output_format}")
