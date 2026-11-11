"""CSV serialization of selection results."""

import csv
import json
import os


class CsvResultWriter:
    """Write one CSV row per match, preserving records without matches."""

    def write(self, result, destination, *, sdf_file=None):
        writer = csv.writer(destination, lineterminator="\n")
        writer.writerow([
            "csv_line", "sdf_file", "record_index", "molecule_name", "selector",
            "smarts", "atom_index_base", "match_index", "atom_indices", "site_atom_indices",
        ])
        if "rows" in result:
            groups = (
                (row["csv_line"], row["sdf_file"], item)
                for row in result["rows"] for item in row["results"]
            )
        else:
            groups = [(None, str(sdf_file), result)]
        for line, path, item in groups:
            path = os.path.relpath(path)
            for molecule in item["molecules"]:
                matches = molecule["matches"]
                common = [
                    line, path, molecule["record_index"], molecule["name"],
                    item["selector"]["name"], item["selector"]["smarts"], item["atom_index_base"],
                ]
                if not matches:
                    writer.writerow(common + ["", "", ""])
                for index, match in enumerate(matches):
                    writer.writerow(common + [
                        index,
                        json.dumps(match["atom_indices"], separators=(",", ":")),
                        json.dumps(match["site_atom_indices"], separators=(",", ":")),
                    ])
