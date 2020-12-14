FROM registry-gitlab.pasteur.fr/gem/macsyfinder/run_dep

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

RUN apt install -y make graphviz libgraphviz-dev python3-dev python3-pygraphviz inkscape &&\
    python3 -m pip install sphinx==3.3.1 sphinx-rtd-theme==0.5.0 sphinx_autodoc_typehints==1.7.0 sphinxcontrib-svg2pdfconverter==1.1.0
RUN apt clean -y

CMD ["/bin/bash"]