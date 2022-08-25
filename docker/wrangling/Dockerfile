FROM bitnami/python:3.9.5-debian-10-r48

RUN apt update && apt -y install git python3-pip python3-wheel && pip3 install wheel
RUN apt -y --install-recommends install dirmngr software-properties-common apt-transport-https \ 
	&& apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' \
	&& add-apt-repository 'deb http://cloud.r-project.org/bin/linux/debian buster-cran40/' \
	&& apt update \
	&& apt -y install libssl-dev libcurl4-openssl-dev libxml2-dev \
	&& apt -y install r-base

WORKDIR /src
COPY . .

# install dhart-utils (dependencies) 
RUN pip3 install -r dependencies/requirements.txt
# install R dependencies
RUN Rscript dependencies/R-dependencies.r

ENTRYPOINT [ "python3", "wrangling.py" ]
