FROM ubuntu:zesty

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

RUN apt-get update -y &&\
    apt-get install -y --no-install-recommends python ncbi-blast+ hmmer &&\
    apt-get clean -y

CMD ["/bin/bash"]
