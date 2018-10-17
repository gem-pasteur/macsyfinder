FROM ubuntu:bionic

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y &&\
    apt-get install -y --no-install-recommends python3.6 ncbi-blast+ hmmer &&\
    apt-get clean -y

ENV DEBIAN_FRONTEND teletype
ENV PYTHONIOENCODING UTF-8

CMD ["/bin/bash"]
