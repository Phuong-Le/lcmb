# modified from https://github.com/AlexandrovLab/SigProfilerMatrixGenerator/blob/master/Dockerfile
# Use an official Python runtime as a parent image
# Start with a base Ubuntu image and install Python
FROM ubuntu:22.04

# Avoid prompts from apt
ARG DEBIAN_FRONTEND=noninteractive

# Install Python and other dependencies, and apply updates
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y python3-pip python3-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*


# Install SigProfilerMatrixGenerator from PyPI
RUN pip install SigProfilerMatrixGenerator==1.2.25 \
    && SigProfilerMatrixGenerator install GRCh37 \
    && SigProfilerMatrixGenerator install GRCh38


# Test that SigProfilerMatrixGenerator is installed and can be called 
RUN SigProfilerMatrixGenerator matrix_generator --help