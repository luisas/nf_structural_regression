FROM debian:jessie

RUN apt-get update && \
    apt-get -y install gcc g++ wget --force-yes

RUN wget https://zhanggroup.org/TM-score/TMscore.cpp

RUN g++ -static -O3 -ffast-math -lm -o TMscore TMscore.cpp

RUN mv TMscore /bin/
