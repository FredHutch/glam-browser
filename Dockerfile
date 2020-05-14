FROM python:3.8.2-slim
ADD requirements.txt /home/dash/
RUN apt-get update
RUN apt update && \
	apt install -y hdf5-tools libhdf5-dev libhdf5-serial-dev nginx supervisor redis-server && \
	pip3 install -r /home/dash/requirements.txt && \
	HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/ pip3 install tables
RUN useradd -u 5555 -m -d /home/dash -c "dash user" dash
ADD system/. /home/dash/system/
ADD app.py /home/dash/
ADD helpers/ /home/dash/helpers/
RUN chown -R dash:dash /home/dash 
WORKDIR /home/dash
EXPOSE 8050
ENV DATA_FOLDER=/share
CMD redis-server & python3 /home/dash/app.py
