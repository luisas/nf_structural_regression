FROM debian:jessie

RUN apt-get update && \
    apt-get -y install gcc g++ gfortran && \
    apt-get install -y --no-install-recommends \
      ed \
      less \
      vim-tiny \
      wget \
      git \
      python \
      build-essential \
      cmake \
      curl \
      libargtable2-0 \
      ca-certificates \
      libpng12-0 \
      python-biopython \
      python-numpy \
      ruby \
      python-setuptools \
      default-jdk \
      valgrind \
      libpng-dev \
      r-base \
      r-recommended \
      autoconf automake \
    	make \
    	autoconf automake \
    	tar \
    	bzip2 \
    	libssl-dev \
    	flex \
    	zlib1g-dev \
    	libbz2-dev \
    	libiomp-dev


#install R packages
RUN R -e "install.packages('seqinr', repos='http://R-Forge.R-project.org')"
RUN R -e "install.packages('ape', repos = 'http://ape-package.ird.fr/')"

# install argtable 2 (needed for clustal Omega)
RUN wget http://prdownloads.sourceforge.net/argtable/argtable2-13.tar.gz && \
    tar -zxf argtable2-13.tar.gz && \
    cd argtable2-13 && \
    ./configure && \
    make && \
    make install && \
    rm /argtable2-13.tar.gz


# install MAFFT
RUN wget http://mafft.cbrc.jp/alignment/software/mafft-7.397-with-extensions-src.tgz --no-check-certificate && \
    tar xfvz mafft-7.397-with-extensions-src.tgz && \
    cd mafft-7.397-with-extensions/core/ && \
    sed -i "s/PREFIX = \/usr\/local/PREFIX = \/mafft/g" Makefile && \
    sed -i "s/BINDIR = \$(PREFIX)\/bin/BINDIR = \/mafft\/bin/g" Makefile && \
    make clean && \
    make && \
    make install && \
    wget http://mafft.cbrc.jp/alignment/software/newick2mafft.rb --no-check-certificate && \
    chmod +x newick2mafft.rb && \
    cp newick2mafft.rb /usr/bin/ &&\
    export "PATH=$PATH:/mafft/bin" && \
    export MAFFT_BINARIES="" && \
    cp /mafft/bin/* /usr/bin/. && \
    mv /mafft-7.397-with-extensions /mafft && \
    rm /mafft-7.397-with-extensions-src.tgz

# Install PROBCONS
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    probcons


# Install BLAST+
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/ncbi-blast-2.10.1+-x64-linux.tar.gz && \
    tar -xvf ncbi-blast-2.10.1+-x64-linux.tar.gz && \
    rm ncbi-blast-2.10.1+-x64-linux.tar.gz && \
    cp ncbi-blast-2.10.1+/bin/* /usr/bin/.

# Install PDB DDBB
RUN mkdir -p /database/pdb && \
    cd /database/pdb && \
    wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt && \
    makeblastdb -in pdb_seqres.txt -input_type fasta -dbtype prot

# install FAMSA
RUN git clone https://github.com/refresh-bio/FAMSA.git && \
  cd FAMSA && \
  git checkout e4bffba2be35d48f6c5bf6351a4707c7a94cc628 && \
  sed -i 's/all: famsa famsa-gpu/all: famsa/' makefile && \
  sed -i 's/NO_AVX = false/NO_AVX = true/' makefile && \
  sed -i 's/NO_AVX2 = false/NO_AVX2 = true/' makefile && \
  sed -i 's/NO_GPU = false/NO_GPU = true/' makefile && \
  make && \
  cp famsa /bin/.


# Install TCOFFEE
RUN git clone https://github.com/cbcrg/tcoffee.git tcoffee && \
    cd tcoffee && \
    git checkout 8908d9fc0fd99749be3dcc279fa640b42fda278a && \
    cd t_coffee/src && \
    make t_coffee && \
    cp t_coffee /bin/.

# Install CLUSTAL OMEGA
RUN wget http://www.clustal.org/omega/clustal-omega-1.2.4.tar.gz && \
    tar -zxf clustal-omega-1.2.4.tar.gz && \
    cd clustal-omega-1.2.4 && \
    sed -i '1157s/1024/100/' src/clustal/muscle_tree.c && \
    ./configure && \
    make && \
    make install && \
    rm /clustal-omega-1.2.4.tar.gz
#
# Install SAP
#
RUN git clone https://github.com/mathbio-nimr-mrc-ac-uk/SAP.git && \
    cd SAP/ && git checkout 1d956d950d4e6e848878b4732f7e42e2b414c34c && \
    autoreconf -fvi && ./configure && make && make install && \
    cd ../ && rm -rf SAP/
#
# Install TMalign and TMscore
#
RUN cd / && wget --no-check-certificate https://zhanglab.ccmb.med.umich.edu/TM-align/TMalign.f && \
        gfortran -static -O3 -ffast-math -lm -o TMalign TMalign.f  && \
        mv TMalign /bin/ && rm TMalign.f && \
        wget --no-check-certificate https://zhanglab.ccmb.med.umich.edu/TM-score/TMscore.f && \
        gfortran -static -O3 -ffast-math -lm -o TMscore TMscore.f && \
        mv TMscore /bin/ && rm TMscore.f



# Download scripts
RUN wget https://raw.githubusercontent.com/edgano/datasets-test/scripts/replace_U.pl && \
   chmod 777 replace_U.pl && \
   mv replace_U.pl /bin/.


ENV HOME '/'

ENV FASTA_4_MAFFT '/bin/fasta35'

ENV HOME_4_TCOFFEE '/tmp/.tcoffee/'
ENV CACHE_4_TCOFFEE '/tmp/.tcoffee/cache'
ENV LOCKDIR_4_TCOFFEE '/tmp/.tcoffee/lock'
ENV TMP_4_TCOFFEE '/tmp/.tcoffee/tmp'
ENV TRASH_4_TCOFFEE '/tmp/.tcoffee/trash'

#unset MAFFT_BINARIES
ENV MAFFT_BINARIES ''

#configure perl5 to BLAST templates
ENV PERL5LIB '/usr/share/perl5:/tcoffee/perl/lib/perl5'
