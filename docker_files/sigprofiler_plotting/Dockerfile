# modified from https://github.com/AlexandrovLab/SigProfilerPlotting/blob/master/Dockerfile
# Start with a base Ubuntu image and install Python
FROM ubuntu:22.04

# Avoid prompts from apt
ARG DEBIAN_FRONTEND=noninteractive

# Install Python and other dependencies, and apply updates
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y python3-pip python3-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*


# Set the working directory in the container
WORKDIR /usr/src/app

# Install SigProfilerPlotting using pip
RUN pip3 install pandas==1.5.3 numpy==1.26.1 matplotlib==3.8.2 pillow==10.1.0

RUN pip3 install SigProfilerPlotting==1.3.24
RUN pip3 install pandas==1.5.3 numpy==1.26.1 matplotlib==3.8.2 pillow==10.1.0

RUN pip3 install SigProfilerPlotting==1.3.24

# Test that SigProfilerMatrixGenerator is installed and can be called 
RUN SigProfilerPlotting plotSBS --help