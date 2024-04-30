FROM --platform=linux/amd64 ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ 'Europe/Zurich'

RUN apt-get update

# install R 
RUN apt-get -y install r-base libssl-dev libsasl2-dev


# install required R libraries

## Install labelSeg package using remotes
RUN R -e "install.packages('remotes', repos = 'https://cloud.r-project.org/')"
RUN R -e "remotes::install_github('baudisgroup/labelSeg')"

## Other packages
RUN R -q -e "install.packages(c('R.utils','MASS'),repos = 'https://cloud.r-project.org/')"

WORKDIR /