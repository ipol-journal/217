# use one of the images from this repository: https://github.com/centreborelli/ipol-docker-images/
FROM registry.ipol.im/ipol:v1-py3.9

# copy the code to $bin
ENV bin /workdir/bin
RUN mkdir -p $bin
WORKDIR $bin
COPY . .


WORKDIR $bin/MissStereo/src
RUN cmake CMakeLists.txt
RUN make
WORKDIR $bin/MissStereo/src/bin
RUN cp homography orsa rectify sift size showRect stereoAC selfSimilar subPixel medianFill $bin

WORKDIR $bin/rdpof1D-2.0
RUN make
RUN cp rdpof1D.exe $bin

WORKDIR $bin/rof1D-2.4
RUN make
RUN cp rof1D.exe $bin

RUN apt-get update && apt-get install -y libgsl-dev
WORKDIR $bin/lucas_kanade-1.05
RUN make exe with_openmp=yes
RUN cp lucas_kanade.exe $bin

WORKDIR $bin/imscript
RUN make
WORKDIR $bin/imscript/bin
RUN cp qeasy qauto $bin


# the execution will happen in the folder /workdir/exec
# it will be created by IPOL

# some QoL tweaks
ENV PYTHONDONTWRITEBYTECODE 1
ENV PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION python
ENV PATH $bin:$PATH

# $HOME is writable by the user `ipol`, but 
ENV HOME /home/ipol
# chmod 777 so that any user can use the HOME, in case the docker is run with -u 1001:1001
RUN groupadd -g 1000 ipol && useradd -m -u 1000 -g 1000 ipol -d $HOME && chmod -R 777 $HOME
RUN pip install imageio
USER ipol