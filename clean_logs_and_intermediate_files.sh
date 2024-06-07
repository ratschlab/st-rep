find . |grep -E ".snakemake/|.ipynb_checkpoints|out/tiles|hs_err_pid|.tree|_graph.bin|_graph.weights|.bin|/core"|xargs rm -rf
rm -rf `find -type d -name .ipynb_checkpoints`