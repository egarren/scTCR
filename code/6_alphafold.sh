

module load alphafold/2.1.1 gcc/6.2.0 cuda/11.2

alphafold.py --fasta_paths=/alphafold/tcr/e.fasta \
--is_prokaryote_list=false \
--max_template_date=2022-01-01 \
--db_preset=full_dbs \
--model_preset=multimer \
--output_dir=/alphafold/output/ \
--data_dir=/n/shared_db/alphafold/


  