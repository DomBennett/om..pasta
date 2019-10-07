# Example host distro
FROM ubuntu:latest

# set up PASTA
RUN apt-get update && apt-get -y upgrade

RUN apt-get install -y python3 python3-setuptools default-jre git libgomp1

RUN ln -s $(which python3) /usr/local/bin/python

RUN cd /usr/local/bin

RUN git clone --depth=10 https://github.com/smirarab/pasta.git

RUN git clone --depth=10 https://github.com/smirarab/sate-tools-linux.git

RUN cd pasta && python3 setup.py develop

RUN rm /usr/local/bin/run_pasta_gui.py

# What does this do?
ENV CONTRALIGN_DIR /usr/local/bin/sate-tools-linux

# outsider *requires* working_dir
RUN mkdir /working_dir
WORKDIR /working_dir
