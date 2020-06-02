#!/bin/bash

set -e

DATA_FOLDER=data gunicorn --workers 4 --worker-class gevent --bind 0.0.0.0:8050 app:server
