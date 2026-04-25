FROM rocker/shiny:4.3.2

# System dependencies for sf, rgdal, leaflet spatial packages
RUN apt-get update && apt-get install -y --no-install-recommends \
  libgdal-dev \
  libgeos-dev \
  libproj-dev \
  libudunits2-dev \
  libssl-dev \
  libcurl4-openssl-dev \
  libxml2-dev \
  && rm -rf /var/lib/apt/lists/*

# Install R package dependencies
RUN Rscript -e "\
  install.packages(c( \
    'golem', 'shiny', 'shinydashboard', \
    'dplyr', 'tidyr', 'readr', 'stringr', \
    'leaflet', 'leaflet.extras', \
    'plotly', 'DT', \
    'sf', 'rgbif', \
    'config', 'rlang', \
    'remotes' \
  ), repos = 'https://cloud.r-project.org', Ncpus = 4) \
"

# Copy the package and install it
WORKDIR /srv/shiny-server/edna-viz
COPY . .

RUN Rscript -e "remotes::install_local('.', dependencies = FALSE)"

# Shiny server config
RUN echo 'run_as shiny;\n\
server {\n\
  listen 3838;\n\
  location / {\n\
    app_dir /srv/shiny-server/edna-viz;\n\
    log_dir /var/log/shiny-server;\n\
  }\n\
}' > /etc/shiny-server/shiny-server.conf

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
