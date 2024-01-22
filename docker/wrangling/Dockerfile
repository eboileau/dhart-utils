FROM bitnami/python:3.10-debian-11

RUN apt update && apt -y install git python3-pip python3-wheel
RUN pip install --upgrade pip wheel setuptools
RUN apt -y --install-recommends install dirmngr software-properties-common apt-transport-https \ 
	&& apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' \
	&& add-apt-repository 'deb http://cloud.r-project.org/bin/linux/debian bullseye-cran40/' \
	&& apt update \
	&& apt -y install libssl-dev libcurl4-openssl-dev libxml2-dev \
	&& apt -y install r-base

WORKDIR /dhart-utils
COPY . .

# Install dhart-utils and dependencies 
RUN pip install .
RUN Rscript dependencies/R-dependencies.r

ENTRYPOINT [ "/opt/bitnami/python/bin/wrangling" ]
