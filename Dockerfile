FROM debian:jessie

# Setup build environment
RUN apt-get clean && apt-get update && apt-get install -y --fix-missing \
  bzip2 \
  fftw3-dev \
  gfortran \
  git \
  liblapack-dev \
  libmpich2-dev \
  make \
  mpich2 \
  wget \
&& rm -rf /var/lib/apt/lists/*

# Switch to user-space
ENV HOME /home/nek

# Pull in gslib
WORKDIR /home/nek
RUN git clone -b async https://github.com/maxhutch/gslib.git
WORKDIR /home/nek/gslib
RUN make

# More stuff for testing
WORKDIR /home/nek
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -p $HOME/miniconda
ENV PATH $HOME/miniconda/bin:$PATH
RUN hash -r
RUN conda config --set always_yes yes --set changeps1 no
RUN conda update -q conda
RUN conda info -a
RUN conda install numpy pytest chest cloudpickle
RUN pip install graphviz mapcombine
RUN pip install git+https://github.com/maxhutch/dask --upgrade
RUN pip install git+https://github.com/maxhutch/nekpy --upgrade

# The tests and tools
RUN git clone https://github.com/maxhutch/nek-tools.git
RUN git clone https://github.com/maxhutch/nek-analyze.git 

# Pull in the sources
WORKDIR /home/nek/NekBox
ADD . . 
WORKDIR /home/nek/NekBox/tests

# Run some tests
ENTRYPOINT ["py.test", "-v", "-s"] 
