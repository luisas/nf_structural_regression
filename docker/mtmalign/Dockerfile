FROM debian:jessie

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

# install mtmalign
RUN wget https://yanglab.nankai.edu.cn/mTM-align/version/mTM-align.tar.bz2 && \
  bzip2 -d mTM-align.tar.bz2 && \
  tar -xvf mTM-align.tar && \ 
  cd mTM-align/src &&\
  cp mTM-align /bin/.
