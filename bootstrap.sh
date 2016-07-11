#!/usr/bin/env bash

# Add repos
add-apt-repository ppa:ubuntu-toolchain-r/test
add-apt-repository ppa:george-edison55/cmake-3.x
apt-get update

# Essentials
apt-get install -y build-essential
apt-get install -y git

# Install newest gcc and g++
apt-get install -y gcc-6 g++-6

# Install CMake
apt-get install -y cmake

# Upgrade all
apt-get -y upgrade
