FROM ubuntu:18.04 as base

ENV TZ=Asia/Dubai
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


RUN apt-get update && \
    apt-get -y install gcc g++ gfortran && \
    apt-get install -y --no-install-recommends \
      ed \
      less \
      vim-tiny \
      vim \
      wget \
      git \
      build-essential \
      cmake \
      curl \
      libargtable2-0 \
      ca-certificates \
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


RUN apt-get install -y make build-essential libssl-dev zlib1g-dev \
libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev \
libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev \
libgdbm-dev libnss3-dev libedit-dev libc6-dev

RUN wget https://www.python.org/ftp/python/3.6.15/Python-3.6.15.tgz && \
   tar -xzf Python-3.6.15.tgz &&\
    cd Python-3.6.15 && \
    ./configure --enable-optimizations  -with-lto  --with-pydebug && \
    make -j 8 && \
    make altinstall


# Install TCOFFEE
RUN git clone https://github.com/cbcrg/tcoffee.git tcoffee && \
    cd tcoffee && \
    git checkout ec0be9b62914c24d8c405d3ac7270ff340eab5cd && \
    cd t_coffee/src && \
    make t_coffee && \
    cp t_coffee /bin/.

RUN pip3 install itertools pandas numpy biopython 


#
# Install TMalign and TMscore
#
RUN cd / && wget --no-check-certificate https://zhanglab.ccmb.med.umich.edu/TM-align/TMalign.f && \
        gfortran -static -O3 -ffast-math -lm -o TMalign TMalign.f  && \
        mv TMalign /bin/ && rm TMalign.f && \
        wget --no-check-certificate https://zhanglab.ccmb.med.umich.edu/TM-score/TMscore.f && \
        gfortran -static -O3 -ffast-math -lm -o TMscore TMscore.f && \
        mv TMscore /bin/ && rm TMscore.f
