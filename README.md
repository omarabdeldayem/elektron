# Elektron
Dependency-free, header-only template library for robotics and embedded systems.

## Features
- High portability (as long as you have a valid C++11 compiler)
- Makes strict use of the stack for systems where using `new`/`delete` is not possible. To use the heap, add `#define ELEKTRO_USE_HEAP` just before you include the library.

## To Use
1. Clone the repo to `$YOUR_DIR`.
2. Add `#include '$YOUR_DIR/src/elektro.hpp` where Elektro is needed.

## For Development on Windows
1. Download Vagrant from: https://www.vagrantup.com/downloads.html
2. Clone this repository onto your host machine.
3. From your terminal, cd into the repository and enter `vagrant up`. This will use the `Vagrantfile` in the repository to start up an ubuntu/trusty64 box. It will then provision the box using the `bootstrap.sh` provisioning script.
4. SSH into your box using `vagrant ss` from Windows 10 Bash shell or using an SSH client (PuTTY), with 'vagrant' as the username and password.
