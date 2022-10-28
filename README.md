# scTCR

This repository contains the code used in our single cell sequencing paper: "Loss of B cell Tolerance is TCR Dependent"

## Installation guide
Install dependencies listed below

## Demo
Model data is included in the `data` directory.

## Instructions
1. Download `data` and `code` directories
2. Set `data` as the working directory
3. Download `gex` and `vdj` data from GSE157649, available here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157649
4. Run each script in numerical order in either Rstudio or Python.  These scripts will generate the figures presented in our manuscript.\
\
NB: Expected total run time is 3-5 days

## System requirements and software

AlphaFold v2.1.1
GIANA v4.1
GLIPH v2
tcrdist3 v0.2.2
DeepTCR v2.1.0


R (v4.1.1) packages

 [1] ggalluvial_0.12.3      harmony_0.1.0          Rcpp_1.0.9            
 [4] sp_1.5-0               SeuratObject_4.1.2     Seurat_4.2.0          
 [7] ggseqlogo_0.1          ggforce_0.4.1          phylotools_0.2.4      
[10] ape_5.6-2              viridis_0.6.2          viridisLite_0.4.1     
[13] VennDiagram_1.7.3      futile.logger_1.4.3    limma_3.46.0          
[16] EnhancedVolcano_1.13.2 plotly_4.10.0          biomaRt_2.46.3        
[19] plyr_1.8.7             forcats_0.5.2          stringr_1.4.1         
[22] purrr_0.3.4            readr_2.1.3            tidyverse_1.3.2       
[25] qgraph_1.9.2           cowplot_1.1.1          ggrepel_0.9.1         
[28] ggpubr_0.4.0.999       immunarch_0.8.0        patchwork_1.1.2       
[31] dtplyr_1.2.2           dplyr_1.0.10           data.table_1.14.2     
[34] vegan_2.6-4            lattice_0.20-45        permute_0.9-7         
[37] xlsx_0.6.5             RColorBrewer_1.1-3     scales_1.2.1          
[40] gridExtra_2.3          ggplot2_3.3.6          dendextend_1.16.0     
[43] corrplot_0.92          pheatmap_1.0.12        tidyr_1.2.1           
[46] tibble_3.1.8     


Python (v3.7.4) packages (TCRdist venv)

absl-py==0.11.0
airr==1.3.1
aniso8601==9.0.1
anndata==0.7.6
anndata2ri==1.0.2
astor==0.8.1
astunparse==1.6.3
attrs==19.3.0
backcall==0.1.0
biopython==1.76
bleach==3.3.0
boto3==1.25.1
botocore==1.28.1
cachetools==4.1.1
CellPhoneDB==3.1.0
certifi==2020.12.5
cffi==1.14.0
changeo==1.0.2
chardet==3.0.4
charset-normalizer==2.0.7
click==7.1.2
colorama==0.4.4
compass-sc==0.9.9.6.2
cplex==20.1.0.0
cryptography==38.0.1
cycler==0.10.0
Cython==0.29.23
decorator==4.4.2
DeepTCR==2.1.27
defusedxml==0.6.0
descartes==1.1.0
dill==0.3.6
distinctipy==1.2.1
docutils==0.17
entrypoints==0.3
fbpca==1.0
feather-format==0.4.1
fisher==0.1.10
fishersapi==0.5
Flask==1.0.4
Flask-RESTful==0.3.9
Flask-Testing==0.7.1
flatbuffers==2.0
gast==0.3.3
geosketch==0.3
get_version==2.1
google-auth==1.23.0
google-auth-oauthlib==0.4.2
google-pasta==0.2.0
grpcio==1.34.0
h5py==2.10.0
hierdiff==0.9
idna==2.10
igraph==0.9.8
importlib-metadata==3.10.0
ipykernel==5.2.1
ipython==7.13.0
ipython-genutils==0.2.0
ipywidgets==7.5.1
itsdangerous==2.0.1
jedi==0.17.0
jeepney==0.8.0
Jinja2==2.11.2
jmespath==0.10.0
joblib==0.17.0
jsonschema==3.2.0
jupyter==1.0.0
jupyter-client==6.1.3
jupyter-console==6.1.0
jupyter-core==4.6.3
keras==2.7.0
Keras-Applications==1.0.8
Keras-Preprocessing==1.1.2
keyring==23.0.1
kiwisolver==1.3.1
legacy-api-wrap==1.2
leidenalg==0.8.8
libclang==12.0.0
llvmlite==0.36.0
logomaker==0.8
loompy==3.0.6
Markdown==3.1.1
MarkupSafe==1.1.1
matplotlib==3.3.3
mistune==0.8.4
mizani==0.7.3
more-itertools==8.4.0
mpmath==1.1.0
natsort==7.0.1
nbconvert==5.6.1
nbformat==5.0.6
networkx==2.4
nose==1.3.7
notebook==6.0.3
numba==0.53.1
numexpr==2.7.1
numpy==1.16.5
numpy-groupies==0.9.13
oauthlib==3.1.0
olga==1.2.4
opt-einsum==3.3.0
packaging==20.9
palettable==3.3.0
palmotif==0.4
pandas==1.1.4
pandocfilters==1.4.2
parasail==1.3.3
parmap==1.6.0
parso==0.7.0
patsy==0.5.3
pexpect==4.8.0
pickleshare==0.7.5
pika==1.3.1
Pillow==8.0.1
pkginfo==1.7.0
plotnine==0.8.0
pluggy==0.13.1
presto==0.6.2
progress==1.6
prometheus-client==0.7.1
prompt-toolkit==3.0.5
protobuf==3.11.2
psutil==5.6.7
ptyprocess==0.6.0
pwseqdist==0.6
py==1.9.0
pyarrow==9.0.0
pyasn1==0.4.8
pyasn1-modules==0.2.8
pycparser==2.20
pyfasta==0.5.2
Pygments==2.8.1
pynndescent==0.5.7
pyparsing==2.4.6
pyrsistent==0.16.0
pysam==0.16.0.1
pytest==5.4.3
python-dateutil==2.8.1
python-igraph==0.9.8
python-libsbml==5.19.0
pytz==2019.3
PyYAML==5.4.1
pyzmq==19.0.0
qtconsole==4.7.3
QtPy==1.9.0
readme-renderer==29.0
requests==2.25.0
requests-oauthlib==1.3.0
requests-toolbelt==0.9.1
rfc3986==1.4.0
rpy2==3.5.5
rsa==4.6
s3transfer==0.6.0
scanpy==1.5.1
scikit-learn==1.0.1
scipy==1.7.3
scvelo==0.2.3
seaborn==0.11.0
SecretStorage==3.3.3
Send2Trash==1.5.0
setuptools-scm==4.1.2
simplegeneric==0.8.1
six==1.13.0
sklearn==0.0
SQLAlchemy==1.3.24
statsmodels==0.13.2
svgwrite==1.4.3
sympy==1.5.1
tables==3.6.1
tcrdist3==0.2.2
tcrmatch==0.0.1
tcrsampler==0.1.9
tensorboard==2.7.0
tensorboard-data-server==0.6.1
tensorboard-plugin-wit==1.8.0
tensorflow==2.7.0
tensorflow-estimator==2.7.0
tensorflow-io-gcs-filesystem==0.22.0
termcolor==1.1.0
terminado==0.8.3
testpath==0.4.4
texttable==1.6.2
threadpoolctl==2.1.0
tornado==6.0.4
tqdm==4.60.0
traitlets==4.3.3
twine==3.4.1
typing-extensions==3.7.4.3
tzlocal==2.1
umap-learn==0.5.3
urllib3==1.26.2
velocyto==0.17.17
virtualenv==16.7.4
wcwidth==0.1.9
webencodings==0.5.1
Werkzeug==0.16.0
widgetsnbextension==3.5.1
wrapt==1.12.1
xlrd==1.2.0
yamlordereddictloader==0.4.0
zipdist==0.1.5
zipp==3.4.1



