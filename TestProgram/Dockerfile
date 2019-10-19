FROM centos/devtoolset-7-toolchain-centos7

USER 0
RUN yum install -y --setopt=tsflags=nodocs make perl openmpi3-devel && \
    yum clean all -y && \
    echo export PATH=$PATH:/usr/lib64/openmpi3/bin/ >> /etc/profile.d/ompi3.sh && \
    sed -i 's/bash/bash -l/g' /usr/bin/container-entrypoint
USER 1001
