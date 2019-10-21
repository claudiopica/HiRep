FROM centos/devtoolset-7-toolchain-centos7 as createcache

USER 0
RUN yum install -y epel-release && \
    yum install -y --setopt=tsflags=nodocs make perl openmpi3-devel ccache && \
    yum clean all -y && \
    echo export PATH=$PATH:/usr/lib64/openmpi3/bin/ >> /etc/profile.d/ompi3.sh && \
    sed -i 's/bash/bash -l/g' /usr/bin/container-entrypoint

RUN yum install -y git && yum clean all -y

RUN mkdir -p /github/workspace &&\ 
    cd /github/workspace &&\
    git clone --depth 1 https://github.com/claudiopica/HiRep.git --branch HiRep-COMP --single-branch

RUN cd  /github/workspace/HiRep/TestProgram &&\
    ../Make/Utils/write_mkflags.pl -f ../Make/MkFlags -ccache -no-mpi -n 2 -r FUND &&\
    ( cd .. && make cleanall ) &&\
    make tests

RUN source /etc/profile &&\
    cd  /github/workspace/HiRep/TestProgram &&\
    ../Make/Utils/write_mkflags.pl -f ../Make/MkFlags -ccache -mpi -n 2 -r FUND &&\
    ( cd .. && make cleanall ) &&\
    make tests

RUN source /etc/profile &&\
    cd  /github/workspace/HiRep/TestProgram &&\
    ../Make/Utils/write_mkflags.pl -f ../Make/MkFlags -ccache -no-mpi -n 3 -r FUND &&\
    ( cd .. && make cleanall ) &&\
    make tests

RUN source /etc/profile &&\
    cd  /github/workspace/HiRep/TestProgram &&\
    ../Make/Utils/write_mkflags.pl -f ../Make/MkFlags -ccache -mpi -n 3 -r FUND &&\
    ( cd .. && make cleanall ) &&\
    make tests

FROM centos/devtoolset-7-toolchain-centos7

USER 0
RUN yum install -y epel-release && \
    yum install -y --setopt=tsflags=nodocs make perl openmpi3-devel ccache && \
    yum clean all -y && \
    echo export PATH=$PATH:/usr/lib64/openmpi3/bin/ >> /etc/profile.d/ompi3.sh && \
    sed -i 's/bash/bash -l/g' /usr/bin/container-entrypoint

COPY --from=createcache /opt/app-root/src/.ccache /opt/app-root/src/.ccache

RUN chown 1001:1001 -R /opt/app-root/src/.ccache

USER 1001

RUN ccache --zero-stats