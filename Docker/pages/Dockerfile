FROM registry-gitlab.pasteur.fr/gem/macsyfinder:run_dep

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

RUN apt-get install -y python-sphinx python-pygraphviz make &&\
    apt-get clean -y

CMD ["/bin/bash"]