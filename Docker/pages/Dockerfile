FROM registry-gitlab.pasteur.fr/gem/macsyfinder/run_dep

MAINTAINER Bertrand Neron <bneron@pasteur.fr>

USER root

RUN apt install -y make graphviz libgraphviz-dev python3-dev python3-pygraphviz inkscape &&\
    python3 -m pip install sphinx sphinx-rtd-theme sphinx_autodoc_typehints sphinxcontrib-svg2pdfconverter
RUN apt clean -y

CMD ["/bin/bash"]