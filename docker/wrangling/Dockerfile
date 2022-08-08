FROM bitnami/python:3.9.5-debian-10-r48







# Replace these variables with fitting ones
ARG destination
ARG geo_identifier

ENV destination=${destination}
ENV geo_identifier=${geo_identifier}
ENV geo_bulk_file=''
ENV data_type='single-cell RNA-Seq'
ENV annotation='Ensembl'






# Update the package repository
RUN apt update

# Install needed packages
RUN pip install rpy2

# Install git to be able to clone the DHART repo in the next step
RUN apt install -y git

# Clone the DHART library
RUN git clone https://github.com/eboileau/dhart-utils.git



ENTRYPOINT [ "python3", "/wrangling.py" ]