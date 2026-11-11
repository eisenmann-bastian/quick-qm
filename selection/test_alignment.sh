rm -r test_out_aligned.csv
rm -r test_out_aligned.sdf

python3 alignment.py test_out.csv -o test_out_aligned.csv
python3 alignment.py test_out.csv -o test_out_aligned.sdf --format sdf