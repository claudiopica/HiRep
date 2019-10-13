FROM spack/centos7
LABEL maintainer="Claudio Pica <pica@cp3.sdu.dk>" version="1.0"

WORKDIR /root
SHELL ["/bin/bash", "-l", "-c"]

# install compilers and mpi libraries
RUN source /etc/profile.d/spack.sh \
 && spack install gcc \
 && spack load gcc \
 && spack install openmpi

ENTRYPOINT ["/bin/bash", "/opt/spack/share/spack/docker/entrypoint.bash"]
CMD ["docker-shell"]