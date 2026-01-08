find . \
  \( -type d \( -name ".snakemake" -o -name ".ipynb_checkpoints" -o -name "lightning_logs" \) -prune -print0 \) -o \
  \( -type f -path "*/logs/slurm-*" -print0 \) |
pv -0 -l | xargs -0 -P "$(nproc)" rm -rf



