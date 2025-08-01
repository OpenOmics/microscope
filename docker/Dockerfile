# Base image for microscope, uses 
# rocker/shiny:v4.4.3 as base image
# which uses Ubuntu/24.04 (noble)
FROM rocker/shiny:4.4.3

############### ABOUT ###############
# Depedencies of microscope:
# • https://github.com/OpenOmics/microscope
# @TLDR:
# • R already satisfied: 4.4.3
#     • shiny,   installed via apt-get
#     • ggplot2, installed via apt-get
#     • ggpubr,  installed via apt-get
#     • dplyr,   installed via apt-get
#     • tidyr,   installed via apt-get
#     • readr,   installed via apt-get
#     • stringr, installed via apt-get
#     • plotly,  installed via apt-get
#     • bslib,   installed via apt-get
#     • tibble,  installed via apt-get
#     • ape,     installed via apt-get  
#     • vegan,   installed via apt-get
#     • usedist, installed via github
#     • scales,  installed via apt-get
#     • shinyFiles, installed via apt-get
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

# Install CRAN packages
RUN Rscript - <<EOF
# Misc helper functions 
err <- function(...){cat(sprintf(...), sep='\n', file=stderr())}
fatal <- function(...) {err(...); quit(status = 1)}

# Number of CPUs processes 
# to use for parallelizing
# the install of N packages 
use_ncpus <- max(parallel::detectCores()-2, 2)

# CRAN packages,
# add any missing/required
# CRAN packages to the list
# directly below. This script
# will install them if they
# are not already installed.
# Install via: 
#   install.packages('packageName', repos='http://cran.r-project.org')
cran_packages <- c(
    "ape",
    "bslib",
    "devtools",
    "dplyr",
    "ggplot2",
    "ggpubr",
    "plotly",
    "readr",
    "scales",
    "shiny",
    "shinyFiles",
    "stringr",
    "tibble",
    "tidyr",
    "usedist",
    "vegan"
)

# Get list of missing CRAN packages
missing <- setdiff(
    cran_packages, rownames(installed.packages())
)
# Log what will be installed
print(
    paste0("Installing ", length(missing), " missing packages: ", paste(missing, collapse = ", "))
)
# Install missing CRAN packages 
install.packages(
    missing,
    Ncpus = use_ncpus,
    repos = 'http://cran.r-project.org'
)

# Check for any missing packages
# after installation, some times
# install.packages() continues on
# error and does not stop with a
# non-zero exit code.
# Get list of missing CRAN packages
missing <- setdiff(
    cran_packages, rownames(installed.packages())
)
# If there are still missing packages
# exit with a non-zero exit code
if (length(missing) > 0) {
    err("Error: The following CRAN packages are still missing after installation: %s", paste(missing, collapse = ", "))
    fatal("Please check the installation logs for errors and try again!")
} else {
    print("All required CRAN packages are installed.")
}
EOF

# Install microscope from local
# github repository, must copy
# the shiny directory into the
# CWD before running docker build.
COPY shiny/ /srv/shiny-server/microscope
RUN mkdir -p /srv/shiny-server/microscope/app_cache

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
RUN chown shinyuser:shinyuser /srv/shiny-server/microscope/app_cache
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
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/microscope', port=8080, host='0.0.0.0')"]
