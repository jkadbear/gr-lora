FROM ubuntu:20.04

WORKDIR /lib

RUN sed -i 's/[a-z]*.ubuntu.com/mirror.tuna.tsinghua.edu.cn/g' /etc/apt/sources.list

RUN apt update && apt install -y software-properties-common

RUN add-apt-repository ppa:gnuradio/gnuradio-releases-3.8

RUN apt update

# Install required dependencies
RUN apt install -y gnuradio cmake libboost-all-dev swig liborc-dev python3-sphinx doxygen libsdl1.2-dev libgsl-dev libqwt-qt5-dev libqt5opengl5-dev libzmq3-dev gobject-introspection gir1.2-gtk-3.0

# Install gr-lora
workdir /src
COPY . .
run mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr .. && \
    make -j`nproc` && \
    make install && \
    ldconfig 

workdir /src/examples

expose 52002
