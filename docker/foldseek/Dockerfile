FROM ubuntu:16.04

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


RUN wget https://mmseqs.com/foldseek/foldseek-linux-sse41.tar.gz; tar xvzf foldseek-linux-sse41.tar.gz;  cp ./foldseek/bin/foldseek /bin/.



RUN apt-get update && apt-get install -y zlib1g-dev build-essential wget

RUN apt-get update && apt-get -y upgrade

RUN apt-get install -y python3-pip

RUN pip3 install --upgrade "pip < 21.0"

RUN pip3 install pandas biopython



# Install TCOFFEE
RUN git clone https://github.com/luisas/tcoffee.git tcoffee && \
    cd tcoffee && \
    git checkout 090eb75b73e562d6e1ad82a365965be23a12e7b8 && \
    cd t_coffee/src && \
    make t_coffee && \
    cp t_coffee /bin/.