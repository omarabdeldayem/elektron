#!/usr/bin/env bash

# Add repos
add-apt-repository ppa:ubuntu-toolchain-r/test
add-apt-repository ppa:webupd8team/java -y
apt-get update

# Essentials
apt-get install -y build-essential
apt-get install -y git

# Install newest gcc and g++
apt-get install -y gcc-6 g++-6

# Install and setup Java 8
echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | sudo /usr/bin/debconf-set-selections
apt-get install oracle-java8-installer
echo "Setting environment variables for JAVA8"
apt-get install -y oracle-java8-set-default

# install unzip for Gradle installation
apt-get install -y unzip

# Extract Gradle
wget "https://services.gradle.org/distributions/gradle-2.14-all.zip"
unzip gradle-2.14-all.zip

# Configure git
git config --global user.name "Omar Abdeldayem"
git config --global user.email "omar.abdeldayem@outlook.com"

