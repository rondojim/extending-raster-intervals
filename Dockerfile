# Use an official Ubuntu as a parent image
FROM ubuntu:20.04

# Set the working directory in the container
WORKDIR /usr/src/app

# Preconfigure the timezone to avoid interactive prompts
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London

# Install dependencies
RUN apt-get update && apt-get install -y \
    tzdata \
    build-essential \
    cmake \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

# Set the timezone
RUN ln -fs /usr/share/zoneinfo/$TZ /etc/localtime && dpkg-reconfigure --frontend noninteractive tzdata

# Copy the current directory contents into the container at /usr/src/app
COPY . .

# Create build directory
RUN mkdir build
WORKDIR /usr/src/app/build

# Run CMake and build the project
RUN cmake ..
RUN make

# Set the environment variables
ENV BOOST_INCLUDEDIR=/usr/include/boost
ENV BOOST_LIBRARYDIR=/usr/lib

# Default command to run your app
CMD ["./main_executable"]
