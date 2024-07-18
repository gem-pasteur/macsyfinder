FROM python:3.10-slim-bookworm

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update -y &&\
    apt install -y --no-install-recommends hmmer git
RUN apt clean -y

ENV DEBIAN_FRONTEND teletype
ENV PYTHONIOENCODING UTF-8

RUN useradd -m msf
USER msf
WORKDIR /home/msf
ENV PATH="$PATH:/home/msf/.local/bin"

CMD ["/bin/bash"]
