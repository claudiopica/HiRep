FROM alpine:3.16

RUN apk update \
      && apk upgrade \
      && apk add --no-cache gcc g++ perl openmpi-dev git openssh samurai bash libexecinfo-dev libexecinfo

#### ADD DEFAULT USER ####
ARG USER=hirep
ENV USER ${USER}
RUN adduser -D ${USER}
ENV USER_HOME /home/${USER}
RUN chown -R ${USER}:${USER} ${USER_HOME}
USER ${USER}