FROM bitnami/python:3.9.5-debian-10-r48

RUN apt update && apt install -y git python3-pip python3-wheel
RUN pip3 install wheel

WORKDIR /src
COPY . .

# install dhart-utils (dependencies)
RUN pip3 install -r dependencies/requirements.txt

ENTRYPOINT [ "python3", "accession.py" ]
