FROM ubuntu:16.04
MAINTAINER Vlad Saveliev "https://github.com/vladsaveliev"

# Setup a base system
RUN apt-get update && \
    apt-get install -y curl wget make g++ git tar gzip bzip2 build-essential  \
        python2.7-dev python-pip python-virtualenv zlib1g-dev default-jre && \
    apt-get upgrade -y libstdc++6

# TargQC installation
COPY . Fingerprinting
RUN pip install --upgrade setuptools pip && \
    cd Fingerprinting && \
    pip install --upgrade -r requirements.txt && \
    python setup.py develop && \
    cd ..
