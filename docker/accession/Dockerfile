FROM bitnami/python:3.9.5-debian-10-r48

RUN apt update && apt install -y git python3-pip python3-wheel
RUN pip3 install wheel setuptools

WORKDIR /dhart-utils
COPY . .

# Install dhart-utils
# Note: rpy2 in API mode cannot be built without R - not required for accession 
RUN sed -i 's/rpy2==3.5.3//' dependencies/requirements.txt
RUN pip3 install -r dependencies/requirements.txt .

ENTRYPOINT ["/opt/bitnami/python/bin/accession"] 
