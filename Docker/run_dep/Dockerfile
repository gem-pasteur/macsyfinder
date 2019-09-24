FROM ubuntu:cosmic

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y &&\
    apt-get install -y --no-install-recommends hmmer python3.7
RUN apt-get clean -y

RUN apt install -y python3-pip &&\
    python3.7 -m pip install pip


ENV DEBIAN_FRONTEND teletype
ENV PYTHONIOENCODING UTF-8

CMD ["/bin/bash"]
