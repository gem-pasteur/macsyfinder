# How to build a docker image

Each dockerfile is used to building a specific tag.
Like for singularity.

To build an image

    docker build -f <dockerfile.ext> -t macsyfinder:tag .

for instance for local image integron_finder:2.0.rc6

    docker build -f Dockerfile.2.o.rc6 -t macsyfinder:2.0rc6 .
    
to build the macsyfinder:latest in dockerhub 

    docker login --username=<your_docker_login>
    docker build -f Dockerfile.<tag> -t gempasteur/macsyfinderfinder:<tag>> .
    docker push gempasteur/macsyfinderfinder:<tag>
