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
ADD app.py /home/dash/
RUN chown -R dash:dash /home/dash 
WORKDIR /home/dash
EXPOSE 8050
CMD python3 app.py
