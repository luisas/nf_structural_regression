FROM debian:jessie

RUN	apt-get -y update && \
	apt-get -y install gcc g++ gfortran \
	python python-pip \
	wget curl \
	git \
	make \
	autoconf automake \
	tar \
	bzip2 \
	libssl-dev \
	flex \
	zlib1g-dev \
	libbz2-dev \
	libiomp-dev \
	dssp


#
# Install tcoffee
#
RUN	cd / && git clone https://github.com/cbcrg/tcoffee.git && cd tcoffee && \
	git checkout b3e1bb7b7085643bc05c78b8255b74d149263453 && \
	cd t_coffee/src && make t_coffee && \
	mv * /bin


ENV HOME_4_TCOFFEE '/tmp/.tcoffee/'
ENV CACHE_4_TCOFFEE '/tmp/.tcoffee/cache'
ENV LOCKDIR_4_TCOFFEE '/tmp/.tcoffee/lock'
ENV TMP_4_TCOFFEE '/tmp/.tcoffee/tmp'
ENV TRASH_4_TCOFFEE '/tmp/.tcoffee/trash'
ENV EMAIL_4_TCOFFEE="test@gmail.com"
ENV PERL5LIB="/tcoffee/perl/lib/perl5"

#
# Install TMalign and TMscore
#
RUN cd / && wget --no-check-certificate https://zhanglab.ccmb.med.umich.edu/TM-align/TMalign.f && \
        gfortran -static -O3 -ffast-math -lm -o TMalign TMalign.f  && \
        mv TMalign /bin/ && rm TMalign.f && \
        wget --no-check-certificate https://zhanglab.ccmb.med.umich.edu/TM-score/TMscore.f && \
        gfortran -static -O3 -ffast-math -lm -o TMscore TMscore.f && \
        mv TMscore /bin/ && rm TMscore.f

# Install blast legacy
#
RUN	wget --no-check-certificate https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz && \
	tar -zxvf ncbi-blast-2.9.0+-x64-linux.tar.gz && \
        mv ncbi-blast-2.9.0+/bin/* /bin/ && \
        rm -rf ncbi-blast-2.9.0+*
