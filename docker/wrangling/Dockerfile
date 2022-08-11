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
RUN apt install -y dirmngr --install-recommends
RUN apt install -y software-properties-common
RUN apt install -y apt-transport-https

RUN apt --fix-broken install
RUN apt autoremove



RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bullseye-cran40/'

RUN apt-get update
RUN apt-get upgrade
RUN apt-get install r-base-dev

#RUN apt-get update
#RUN apt-get install -y r-base

# Install git to be able to clone the DHART repo in the next step
RUN apt install -y git

# Clone the DHART library
RUN git clone https://github.com/eboileau/dhart-utils.git

# Install R dependencies
COPY /dependencies/R-dependencies.r /dependencies/R-dependencies.r 
COPY /dependencies/requirements.txt /dependencies/requirements.txt

RUN Rscript /dependencies/R-dependencies.r
RUN pip install -r /dependencies/requirements.txt

# copy handler script
RUN mkdir /root/bin
COPY wrangling.py /root/bin/wrangling.py

ENTRYPOINT [ "python3", "/root/bin/wrangling.py" ]