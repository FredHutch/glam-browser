#!/bin/bash

set -e

# DATA_FOLDER=data gunicorn --workers 4 --worker-class gevent --bind 0.0.0.0:8050 --timeout 120 app:server
DATA_FOLDER=data python3 app.py
