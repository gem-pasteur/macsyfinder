FROM ubuntu:cosmic

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y &&\
    apt-get install -y --no-install-recommends hmmer python3.6 python3.7
RUN apt-get install -y python3-pip &&\
    apt-get clean -y

ENV DEBIAN_FRONTEND teletype
ENV PYTHONIOENCODING UTF-8

CMD ["/bin/bash"]
