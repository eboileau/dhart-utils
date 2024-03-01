FROM docker.io/bitnami/python:3.10-debian-11

RUN apt update && apt install -y git python3-pip python3-wheel
RUN pip install --upgrade pip wheel setuptools

WORKDIR /dhart-utils
COPY . .

# Install dhart-utils
# Note: rpy2 in API mode cannot be built without R - not required for accession 
RUN sed -i 's/rpy2//' setup.cfg
RUN pip install .

ENTRYPOINT ["/opt/bitnami/python/bin/accession"] 
