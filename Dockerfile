FROM python:3.8.2-slim
RUN apt-get update && \
	apt-get install -y hdf5-tools libhdf5-dev libhdf5-serial-dev build-essential && \
	apt-get install -y python3-numpy python3-scipy python3-pandas
ADD requirements.txt /home/dash/
RUN pip3 install -r /home/dash/requirements.txt && \
	pip3 install scikit-bio && \
	HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/ pip3 install tables
RUN useradd -u 5555 -m -d /home/dash -c "dash user" dash
ADD app.py /home/dash/
ADD redis.conf /home/dash/
ADD helpers/ /home/dash/helpers/
ADD share/ /share/
RUN chown -R dash:dash /home/dash 
WORKDIR /home/dash
EXPOSE 8050
ENV DATA_FOLDER=/share
CMD gunicorn --timeout 120 --workers 4 --worker-class gevent --bind 0.0.0.0:8050 app:server
