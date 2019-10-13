FROM spack/centos7
LABEL maintainer="Claudio Pica <pica@cp3.sdu.dk>" version="1.0"

WORKDIR /root
SHELL ["/bin/bash", "-l", "-c"]

# install compilers and mpi libraries
RUN spack install gcc

RUN source /etc/profile.d/spack.sh \
 && spack load gcc \
 && spack compiler add \
 && spack install openmpi

RUN source /etc/profile.d/spack.sh \
 && spack load gcc \
 && spack compiler add \
 && spack install mpich

ENTRYPOINT ["/bin/bash", "/opt/spack/share/spack/docker/entrypoint.bash"]
CMD ["docker-shell"]