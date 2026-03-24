import os
base_dir = "" #NOTE: insert the base directory here
ids = [] #NOTE: Fill in the list of ids here or put them in a config.json file
ct_colnames = ["cell_type_level_1", "cell_type_level_2", "cell_type_level_3"] #NOTE: change here if the cell type columns are not these values
rule all:
    input:
        expand(os.path.join("{base_dir}", "Processed_scRNA/data/magma/{id}/{ct_colname}", "place_holder.txt"), base_dir=base_dir, id=ids, ct_colname=ct_colnames)

rule compute_sumstat:
    input:
        os.path.join("{base_dir}", "Preprocessing_scRNA/data/{id}", "{id}.h5ad")
    output:
        os.path.join("{base_dir}", "Processed_scRNA/data/magma/{id}/{ct_colname}", "place_holder.txt")
    params:
        outdir = os.path.join("{base_dir}", "Processed_scRNA/data/magma/{id}"),
        ct_colname = "{ct_colname}"

    shell:
        """
        touch {output};
        python compute_sumstat_magma.py --h5ad {input} --outdir {params.outdir} --ct_colname {params.ct_colname}
        """