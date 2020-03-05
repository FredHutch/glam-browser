# TODO use r-shiny-server for Safari compatibility
FROM fredhutch/r-shiny-base:3.6.2
RUN apt-get update
RUN apt update && \
	apt install -y pandoc build-essential python3 python3-pip && \
	apt install -y hdf5-tools libhdf5-dev libhdf5-serial-dev nginx supervisor && \
	pip3 install pandas==0.24.2 numpy && \
	HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/ pip3 install tables
RUN R -e "install.packages('tidyquant', repos = 'http://cran.us.r-project.org')" && \
	R -e "install.packages('reticulate', repos = 'http://cran.us.r-project.org')" && \
	R -e "install.packages('shinydashboard', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')" && \
	R -e "BiocManager::install('ComplexHeatmap')"
RUN useradd -u 5555 -m -d /home/shiny -c "shiny user" shiny
ADD app/. /home/shiny/
ADD system/. /home/shiny/system/
RUN chown -R shiny:shiny /home/shiny 
WORKDIR /home/shiny
ARG DATA_DIR
ENV DATA_DIR ${DATA_DIR:-/home/shiny}
ARG APP_DIR
ENV APP_DIR ${APP_DIR:-/home/shiny}
EXPOSE 7777
CMD /usr/bin/supervisord -c /home/shiny/system/sup.conf
