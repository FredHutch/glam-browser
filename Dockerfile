# TODO use r-shiny-server for Safari compatibility
FROM python:3.8.2-slim
ADD requirements.txt /home/dash/
RUN apt-get update
RUN apt update && \
	apt install -y hdf5-tools libhdf5-dev libhdf5-serial-dev nginx supervisor && \
	pip3 install -r /home/dash/requirements.txt && \
	HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/ pip3 install tables
RUN useradd -u 5555 -m -d /home/dash -c "dash user" dash
ADD system/. /home/dash/system/
RUN chown -R dash:dash /home/dash 
WORKDIR /home/dash
ARG DATA_DIR
ENV DATA_DIR ${DATA_DIR:-/home/dash}
ARG APP_DIR
ENV APP_DIR ${APP_DIR:-/home/dash}
ARG PYTHON_BIN
ENV PYTHON_BIN ${PYTHON_BIN:-/usr/bin/python3.7}
EXPOSE 7777
CMD /usr/bin/supervisord -c /home/dash/system/sup.conf
