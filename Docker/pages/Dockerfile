FROM registry-gitlab.pasteur.fr/gem/macsyfinder/run_dep

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

RUN apt-get install -y make graphviz libgraphviz-dev python3.7-dev && apt-get clean &&\
    python3.7 -m pip install sphinx==2.2.0 sphinx-rtd-theme==0.4.3 pygraphviz==1.5

CMD ["/bin/bash"]