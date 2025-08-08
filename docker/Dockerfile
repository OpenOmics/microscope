# Base image for microscope, uses
# rocker/shiny:v4.4.3 as base image
# which uses Ubuntu/24.04 (noble)
FROM rocker/shiny:4.4.3

############### ABOUT ###############
# Depedencies of microscope:
# • https://github.com/OpenOmics/microscope
# @TLDR:
# • R already satisfied: 4.4.3
#     • shiny,   installed via CRAN
#     • ggplot2, installed via CRAN
#     • ggpubr,  installed via CRAN
#     • dplyr,   installed via CRAN
#     • tidyr,   installed via CRAN
#     • readr,   installed via CRAN
#     • stringr, installed via CRAN
#     • plotly,  installed via CRAN
#     • bslib,   installed via CRAN
#     • tibble,  installed via CRAN
#     • ape,     installed via CRAN
#     • vegan,   installed via CRAN
#     • usedist, installed via CRAN
#     • scales,  installed via CRAN
#     • pak,     installed via CRAN
#     • bs4Dash, installed via CRAN
LABEL maintainer=kuhnsa@nih.gov

############### INIT ################
# Create Container filesystem specific
# working directory and opt directories
# to avoid collisions with the host's
# filesystem, i.e. /opt and /data
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2

# Set time zone to US east coast
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
        && echo $TZ > /etc/timezone

############### SETUP ################
# This section installs system packages
# required for your project. If you need
# extra system packages add them here.
# Get system library dependencies of
# R packages by running within R:
# pak::pkg_sysreqs(c("devtools","packageX","packageY","packageZ"), sysreqs_platform = "ubuntu-24.04")
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        bzip2 \
        ca-certificates \
        cmake \
        curl \
        gawk \
        git \
        gzip \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libgit2-dev \
        libharfbuzz-dev \
        libicu-dev \
        libjpeg-dev \
        libpng-dev \
        libssl-dev \
        libtiff-dev \
        libx11-dev \
        libxml2-dev \
        locales \
        make \
        pandoc \
        python3 \
        python3-pip \
        unzip \
        wget \
        zlib1g-dev \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Make python3 the default interpreter
# and install Python Packages
RUN ln -sf /usr/bin/python3 /usr/bin/python

############### INSTALL ################
# Install GNU which, increases compatiability
# with which distros included with Debian/Rocky
# linux installations, this version of which
# includes extra options that the ubuntu version
# does not have. The gnu verison is backwards
# compatiable with the ubunutu version.
RUN mkdir -p /opt2/ \
        && wget --no-check-certificate \
            -O /opt2/which.tar.gz \
            https://ftp.gnu.org/gnu/which/which-2.21.tar.gz \
        && tar -xvf /opt2/which.tar.gz -C /opt2/ \
        && rm /opt2/which.tar.gz \
        && cd /opt2/which-2.21/ \
        && ./configure --prefix=/opt2/which-2.21 \
        && make \
        && make install

# Install microscope shiny application via pak,
# and create app_cache directory for shiny runtime
RUN git clone https://github.com/OpenOmics/microscope.git /srv/shiny-server/microscope \
    && cd /srv/shiny-server/microscope \
    && mkdir -p /srv/shiny-server/microscope/shiny/app_cache \
    && Rscript -e "use_ncpus <- max(parallel::detectCores()-2, 2); \
        if (!requireNamespace('pak', quietly = TRUE)) install.packages('pak', repos = 'https://cloud.r-project.org/', Ncpus = use_ncpus)" \
    && Rscript -e "pak::local_install_dev_deps('/srv/shiny-server/microscope', dependencies = TRUE)"

################ POST #################
# Add Dockerfile and export important
# environment variables
ADD Dockerfile /opt2/microscope.dockerfile
RUN chmod -R a+rX /opt2 /srv/shiny-server
ENV PATH="/opt2/which-2.21/bin:${PATH}"
ENV TMPDIR="/tmp"
WORKDIR /data2
# Create a non-root user, shinyuser,
# to run the Shiny application
RUN useradd -ms /bin/bash shinyuser
RUN chown shinyuser:shinyuser /srv/shiny-server/microscope/shiny/app_cache
USER shinyuser
# Liston on port for incoming connections
EXPOSE 8080
# Add healthcheck to ensure container is
# running and ready to accept connections,
# if not exit with non-zero exit code. This
# is useful for orchestration tools that
# require need a healthcheck to determine
# if the container is ready to accept
# connections. If the healthcheck fails
# after 3 retries the container will be
# automatically restarted.
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
        CMD curl -f http://localhost:8080 \
        || (echo 'Error: Health check failed, unable to reach Shiny app on port 8080!' >&2 && exit 1)
# Run the microscope shiny application
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/microscope/shiny', port=8080, host='0.0.0.0')"]
