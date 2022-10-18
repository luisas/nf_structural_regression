FROM debian:11
LABEL authors="Athanasios Baltzis" \
    title="MMseqs2 docker image" \
    Version="v0.1" \
    description="Docker image containing all software requirements to download the databases and parameteres for Colabfold"

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    wget \
    aria2 \
    procps \
    tar \
    rsync \
    git \
    cmake \
    build-essential \
    && rm -rf /var/lib/apt/lists/*
RUN git clone https://github.com/soedinglab/MMseqs2.git \
&& cd MMseqs2 \
&& git checkout a81d9e727490eca35bea60a53417da6b2428fddb \
&& mkdir build \
&& cd build \
&& cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. .. \
&& make \
&& make install
ENV PATH="/MMseqs2/build/bin:$PATH"
