FROM centos/devtoolset-7-toolchain-centos7

USER 0
RUN yum install -y --setopt=tsflags=nodocs make perl && yum clean all -y
USER 1001
WORKDIR /github/workspace/TestProgram