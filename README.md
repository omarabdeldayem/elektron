# Elektron

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/15d01da01e8f4d20a130c06d8deb70a3)](https://www.codacy.com/app/omarabdeldayem/elektron?utm_source=github.com&utm_medium=referral&utm_content=omarabdeldayem/elektron&utm_campaign=badger)

Dependency-free, header-only template library, suitable for robotics and embedded systems. By default, elektron won't make any heap allocations (other than in-place heap allocations at the start) which is appropriate for small matrices / if you can't use new/delete. However, it's a good idea to add `#define ELEKTRON_USE_HEAP` just before the elektron include for large matrices or if you're planning on using anything in cv or nn (possibly necessary dependeing on how large your matrices are).

## Structure
- **math:** Includes a matrix class with all the expected operations, useful decompostions and a linear equation solver. Also includes numerical differentiation/integration functions.
- **filters:** Kalman and Extended Kalman filters for signal processing.
- **controllers:** Just a standard PID controller thus far.
- **cv:** Currently only includes an Image class for 3-channel images with colour space conversions, a Kernel struct and convolution functions. This isn't intended to be a rewrite of OpenCV - just the essentials.
- **nn:** Mostly empty - use TinyDNN in the meantime. 

## Requirements
A C++11 compiler is all! 

## To Use
1. Clone the repo to `$YOUR_DIR`.
2. Add `#include '$YOUR_DIR/elektron/elektron.hpp` where Elektron is needed.

## For Development on Windows
1. Download Vagrant from: https://www.vagrantup.com/downloads.html
2. Clone this repository onto your host machine.
3. From your terminal, cd into the repository and enter `vagrant up`. This will use the `Vagrantfile` in the repository to start up an ubuntu/trusty64 box. It will then provision the box using `bootstrap.sh`.
4. SSH into your box using `vagrant ss` from Windows 10 Bash shell or using an SSH client (PuTTY), with 'vagrant' as the username and password.
