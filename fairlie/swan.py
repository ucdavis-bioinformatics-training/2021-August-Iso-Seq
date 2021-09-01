import swan_vis as swan

annot_gtf = '/dfs6/pub/freese/mortazavi_lab/ref/gencode.vM21/gencode.vM21.annotation.gtf'
data_gtf = 'all_talon_observedOnly.gtf'
ab_file = 'all_talon_abundance_filtered.tsv'
talon_db = 'c2c12_bulk.db'
pass_list = 'all_pass_list.csv'
meta = 'metadata.tsv'

## from GTF + ab file

sg = swan.SwanGraph()
sg.add_annotation(annot_gtf)
sg.add_transcriptome(data_gtf)
sg.add_abundance(ab_file)
sg.save_graph('swan')

## from DB

sg = swan.SwanGraph()
sg.add_annotation(annot_gtf)
sg.add_transcriptome(talon_db, pass_list=pass_list)
sg.add_abundance(ab_file)
sg.add_metadata(meta)
sg.save_graph('swan')


obs_col = 'time_point'
obs_conditions = ['0hr', '72hr']

# perform a differential gene expression
# Wald test on the provided metadata column and conditions
test = sg.de_gene_test(obs_col, obs_conditions=obs_conditions)
sg.save_graph('swan')


uns_key = swan.make_uns_key('deg',
                            obs_col=obs_col,
                            obs_conditions=obs_conditions)
test = sg.adata.uns[uns_key]
de_genes = sg.get_de_genes(obs_col, obs_conditions=obs_conditions,
                           q=0.05, log2fc=1)

# perform a differential transcript expression
# Wald test on the provided metadata column and conditions
test = sg.de_transcript_test(obs_col, obs_conditions=obs_conditions)
sg.save_graph('swan')

# det - differential transcript expression
uns_key = swan.make_uns_key('det',
                            obs_col=obs_col,
                            obs_conditions=obs_conditions)
test = sg.adata.uns[uns_key]
# return a table of significantly differentially-expressed genes
# for a given q val + log2fc threshold
de_transcripts = sg.get_de_transcripts(obs_col, obs_conditions=obs_conditions,
                           q=0.05, log2fc=1)

# die
die_table = sg.die_gene_test(obs_col=obs_col,
                             obs_conditions=obs_conditions,
                             verbose=True)
sg.save_graph('swan')

es_df = sg.find_es_genes(verbose=True)
sg.save_graph('swan')


ir_df = sg.find_ir_genes(verbose=True)
sg.save_graph('swan')
