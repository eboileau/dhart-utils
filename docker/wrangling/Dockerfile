FROM bitnami/python:3.9.5-debian-10-r48

# Update the package repository
RUN apt update

# Install required packages
RUN pip install rpy2
RUN pip install scipy
RUN pip install pandas
RUN pip install anndata
RUN pip install mygene

# Install R
RUN apt-get update && apt-get install -y r-base

# Install git to be able to clone the DHART repo in the next step
RUN apt install -y git

# Clone the DHART library
RUN git clone https://github.com/eboileau/dhart-utils.git

# Install R dependencies
COPY ./R-dependencies.r ./R-dependencies.r
RUN Rscript ./R-dependencies.r

# copy handler script
RUN mkdir /root/bin
COPY ./wrangling.py /root/bin/wrangling.py

ENTRYPOINT [ "python3", "/root/bin/wrangling.py" ]