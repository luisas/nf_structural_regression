FROM debian:jessie
MAINTAINER Luisa Santus <luisa.santus95@gmail.com>

RUN apt-get update && \
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
      r-recommended

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
