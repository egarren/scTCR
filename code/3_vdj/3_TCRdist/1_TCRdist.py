#see https://github.com/kmayerb/tcrdist3
pip install tcrdist3


#run
import pandas as pd
from tcrdist.repertoire import TCRrep
import numpy as np
from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory2
from tcrdist.rep_funcs import  compute_n_tally_out_of_memory2
from hierdiff.association_testing import cluster_association_test
from tcrdist.public import TCRpublic
import os, sys
os.chdir("/data")


df = pd.read_csv("dash.csv")

tr = TCRrep(cell_df = df,
            organism = 'mouse',
            chains = ['alpha','beta'],
            db_file = 'alphabeta_gammadelta_db.tsv',
            compute_distances = False)

tr.cpus = 10
tr.compute_distances()

tr.pw_alpha
tr.pw_beta
tr.pw_cdr3_a_aa
tr.pw_cdr3_b_aa

np.savetxt("pw_alpha.csv",tr.pw_alpha,delimiter=",")
np.savetxt("pw_beta.csv",tr.pw_beta,delimiter=",")
np.savetxt("pw_cdr3_a_aa.csv",tr.pw_cdr3_a_aa,delimiter=",")
np.savetxt("pw_cdr3_b_aa.csv",tr.pw_cdr3_b_aa,delimiter=",")

##probability of generation
from tcrdist.pgen import OlgaModel

# Load OLGA model as a python object
olga_beta  = OlgaModel(chain_folder = "mouse_T_beta", recomb_type="VDJ")
import parmap
tr.clone_df['pgen_cdr3_b_aa'] = \
    parmap.map(
        olga_beta.compute_aa_cdr3_pgen, 
        tr.clone_df['cdr3_b_aa'], 
        pm_pbar=True, 
        pm_processes = 10)
# Note one can find neighbors based on paired-chain distances.
from tcrdist.public import _neighbors_fixed_radius
from tcrdist.public import _K_neighbors_fixed_radius
tr.clone_df['neighbors'] = _neighbors_fixed_radius(pwmat = tr.pw_beta + tr.pw_alpha, radius = 50)
tr.clone_df['K_neighbors'] = _K_neighbors_fixed_radius(pwmat = tr.pw_beta + tr.pw_alpha , radius = 50)
tr.clone_df['pgen_cdr3_b_aa_nlog10'] = tr.clone_df['pgen_cdr3_b_aa'].apply(lambda x : -1*np.log10(x))
tr.clone_df['nsubject'] = tr.clone_df['neighbors'].apply(lambda x: len(tr.clone_df['subject'][x].unique()))
# nsubject > 1 implies quasi-publicity)
tr.clone_df['qpublic'] = tr.clone_df['nsubject'].apply(lambda x: x >1 )
tr.clone_df.to_csv("nn_vs_k.csv",index=False)


####compare repetoires
import pandas as pd 
from tcrdist.repertoire import TCRrep
from tcrdist.setup_tests import download_and_extract_zip_file

df = pd.read_csv("dash.csv")
df_AID=df.loc[df['epitope'] == 'AID']
df_m564=df.loc[df['epitope'] == 'm564']

# df_AID = pd.read_csv('AID.csv', sep = ",")        #(1)
tr_AID = TCRrep(cell_df = df_AID,                 #(2)
            organism = 'mouse', 
            chains = ['alpha','beta'], 
            db_file = 'alphabeta_gammadelta_db.tsv',
            compute_distances = False)
df_search = tr_AID.clone_df.copy()                 #(3)

tr = TCRrep(cell_df = df_m564,                           #(4)
            organism = 'mouse', 
            chains = ['alpha','beta'], 
            db_file = 'alphabeta_gammadelta_db.tsv')

tr.compute_rect_distances(df = tr.clone_df, df2 = df_search, store =True) #(5)


np.savetxt("rw_beta.csv",tr.rw_beta,delimiter=",")
np.savetxt("rw_cdr3_b_aa.csv",tr.rw_cdr3_b_aa,delimiter=",")


# see how many are tcrdistance less than 50
import numpy as np
np.sum(tr.rw_beta < 50, axis = 1)

###tree diagram
from tcrdist.tree import TCRtree
tcrtree = TCRtree(tcrrep = tr, html_name = 'dash.mouse.ab.tree.html')
tcrtree.build_tree()
assert os.path.isfile('dash.mouse.ab.tree.html')

##hclust
from tcrdist.rep_diff import hcluster_diff
tr.hcluster_df, tr.Z =\
    hcluster_diff(clone_df = tr.clone_df,
                  pwmat    = tr.pw_alpha,
                  x_cols = ['epitope'],
                  count_col = 'count')


#Sankey plots
from tcrdist.plotting import plot_pairings, _write_svg
svg_AID  = plot_pairings(tr.clone_df.loc[tr.clone_df.epitope == "AID"], 
            cols = ['v_b_gene', 'j_b_gene','j_a_gene', 'v_a_gene'], 
            count_col='count')
svg_m564 = plot_pairings(tr.clone_df.loc[tr.clone_df.epitope == "m564"], 
        cols = ['v_b_gene', 'j_b_gene', 'j_a_gene', 'v_a_gene'], 
        count_col='count')
_write_svg(svg_AID, name = "AID_gene_usage_plot.svg", dest = ".")
_write_svg(svg_m564, name = "m564_gene_usage_plot.svg", dest = ".")

#fisher test
import fishersapi
fisherresult=fishersapi.fishers_frame(tr.clone_df.loc[tr.clone_df.epitope == "AID"], 
                         col_pairs=[('v_b_gene', 'j_b_gene'),
                                    ('v_a_gene', 'j_a_gene'),
                                    ('v_a_gene', 'v_b_gene'),
                                    ('j_a_gene', 'j_b_gene')])
fisherresult.to_csv("fisher.csv",index=False)

#public
tp = TCRpublic(
    tcrrep = tr, 
    output_html_name = "quasi_public_clones.html")

public = tp.report()

#fixed radius neighborhood
from tcrdist.rep_diff import neighborhood_diff
# diff testing is pasted on binary comparison, so all epitope not 'PA' are set to 'X'
tr.clone_df['AID'] = ['AID' if x == 'AID' else 'm564' for x in tr.clone_df.epitope]
# Larger Radius
tr.nn_df = neighborhood_diff(clone_df= tr.clone_df, 
    pwmat = tr.pw_beta + tr.pw_alpha, 
    count_col = 'count', 
    x_cols = ['AID'], 
    knn_radius = 150) #adjust as needed
trdistnn=tr.nn_df[['K_neighbors', 'val_0', 'ct_0', 'val_2', 'ct_2', 'RR','OR', 'pvalue', 'FWERp','FDRq']].\
    sort_values(['FDRq']).sort_values(['OR','ct_0'], ascending = False)
trdistnn.to_csv("trdist.csv",index=False)

#hierarchical neighborhood
# diff testing is pasted on binary comparison, so all epitope not 'PA' are set to 'X'
res, Z= hcluster_diff(tr.clone_df, tr.pw_beta, x_cols = ['AID'], count_col = 'count')
res_summary = member_summ(res_df = res, clone_df = tr.clone_df, addl_cols=['epitope'])
res_detailed = pd.concat([res, res_summary], axis = 1)
html = plot_hclust_props(Z,
            title='PA Epitope Example',
            res=res_detailed,
            tooltip_cols=['cdr3_b_aa','v_b_gene', 'j_b_gene','epitope'],
            alpha=0.00001, colors = ['blue','gray'],
            alpha_col='pvalue')
with open('hierdiff_example.html', 'w') as fh:
    fh.write(html)



