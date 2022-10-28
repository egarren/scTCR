#see https://github.com/s175573/GIANA

#install 
pip freeze > requirements.txt
pip install -r requirements.txt --upgrade
pip install numpy --upgrade

#new env
module load gcc/6.2.0 R/4.1.1 python/3.7.4 pandoc/2.1.1 gdal/3.1.4 samtools/1.9 boost/1.75.0  geos/3.8.1 cairo/1.16.0 cmake/3.22.2 miniconda3/4.10.3
unset PYTHONPATH

conda create -n faiss_env -c pytorch python=3.7.4 faiss-cpu
source activate faiss_env


# #install
# pip3 install biopython 
# # conda install -c pytorch faiss-cpu
# conda install -c conda-forge faiss-cpu
# pip3 install -U scikit-learn 
# pip3 install query
# conda install nomkl numpy scipy scikit-learn numexpr
# conda remove mkl


# python3 GIANA4.1.py -f tutorial.txt -N 8 -S 0.1 -t 12
python3 GIANA4.1.py -f tutorial.txt -v -S 10 -N 8
python3 GIANA4.1.py -f tutorial.txt -v -S 10 -N 8 -M
python3 GIANA4.1.py -f AID_tutorial.txt -v -S 10 -N 8
python3 GIANA4.1.py -f m564_tutorial.txt -v -S 10 -N 8


##query
python3 GIANA4.1.py -q m564_tutorial.txt -r AID_tutorial.txt -v -S 10 -N 8
