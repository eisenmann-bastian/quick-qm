#!/bin/sh
set -eu

root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
config_name=${1:-xyz_grid.yaml}
case "$config_name" in
    *[!a-zA-Z0-9_.-]*|.|..|'')
        echo "Invalid example configuration name: $config_name" >&2
        exit 2
        ;;
esac
if [ "$#" -gt 2 ]; then
    echo "Usage: $0 [config-name] [output-directory]" >&2
    exit 2
fi
config="$root/examples/configs/$config_name"
if [ ! -f "$config" ]; then
    echo "Example configuration not found: $config" >&2
    exit 2
fi
if [ "$#" -eq 2 ]; then
    mkdir -p -- "$2"
    output_dir=$(CDPATH= cd -- "$2" && pwd)
else
    output_dir=$(mktemp -d "${TMPDIR:-/tmp}/molecular-aligner-example.XXXXXXXX")
fi

cd "$root/examples"
"${PYTHON:-python3}" "$root/aligner/aligner.py" test_out_aligned.csv \
    --config "$config" --output "$output_dir/poses.sdf" --no-overwrite-output
printf 'Example output: %s\n' "$output_dir"
