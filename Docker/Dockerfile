FROM mambaorg/micromamba:latest

# 1) Initialize micromamba shell (if needed)
RUN micromamba shell init --shell=bash --root-prefix=/opt/conda

# 2) Install spoa and other dependencies
RUN micromamba install -y \
      -p /opt/conda \
      -c conda-forge \
      -c bioconda \
      spoa \
      medaka \
      racon \
      minimap2 \
      blast \
      fastqc \
      multiqc \
      cutadapt \
      nanoplot \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"

# 3) Then install NGSpeciesID via pip
RUN pip install NGSpeciesID

WORKDIR /data
ENTRYPOINT ["bash"]
