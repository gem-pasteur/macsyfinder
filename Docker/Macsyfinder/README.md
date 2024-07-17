# How to build a docker image

Each dockerfile is used to building a specific tag.
Like for singularity.

To build an image

    DOCKER_BUILDKIT=1 docker build -f <dockerfile.ext> -t macsyfinder:tag .

for instance for local image macsyfinder:2.0.rc6

    DOCKER_BUILDKIT=1 docker build -f Dockerfile.2.o.rc6 -t macsyfinder:2.0rc6 .

to build the macsyfinder:latest in dockerhub

    docker login --username=<your_docker_login>
    DOCKER_BUILDKIT=1 docker build -f Dockerfile.<tag> -t gempasteur/macsyfinder:<tag> .
    docker push gempasteur/macsyfinder:<tag>

to build a container to test the local repo

    docker build -f Dockerfile.dep -t msf_test .
    docker run -v ${PWD}:/home/msf  -u $(id -u ${USER}):$(id -g ${USER}) -it msf_test

    python tests/run_tests.py -vv
    ...
