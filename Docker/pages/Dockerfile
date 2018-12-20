FROM registry-gitlab.pasteur.fr/gem/macsyfinder:run_dep

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

RUN apt-get install -y python3-sphinx python3-sphinx-rtd-theme python3-pygraphviz make &&\
    apt-get clean -y

CMD ["/bin/bash"]