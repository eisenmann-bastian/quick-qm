#rm -r ./sample_csv/
#rm -r ./sample_mixed/


# python3 downloader.py sample.csv -o ./sample_csv/
python3 selector.py --csv reference_molecules.csv --format sdf -o reference_combined.sdf

  