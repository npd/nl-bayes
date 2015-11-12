NL-Bayes image denoising
========================

About
-----

* Author    : Marc Lebrun <marc.lebrun.ik@gmail.com>
* Copyright : (C) 2013 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt

This version was modified to allow running just one step of the algorithm, and is compatible with more image formats.

All images are supposed in the range [0, 255].

Overview
--------

This source code provides an implementation of the NL-Bayes image denoising.

Building
--------

To compile, use

    $ mkdir build
    $ cd build
    $ cmake .. [-DCMAKE_CXX_COMPILER=/path/of/c++/compiler -DCMAKE_C_COMPILER=/path/of/c/compiler] [-DCMAKE_BUILD_TYPE=Debug]
    $ make

To rebuild, e.g. when the code is modified, use

    $ cd build
    $ make

Using
-----

To denoise an image `noisy.png` in the range [0, 255] corrupted with a noise of standard deviation `sigma`, use

    $ nl_bayes noisy.png sigma result.tiff

To compute just the first step of the algorithm, use

    $ nl_bayes noisy.png sigma result.tiff -1

To compute just the second step, with `guide.tiff` as guide, use

    $ nl_bayes noisy.png sigma result.tiff -2 guide.tiff

The option `-v` increases verbosity.

To enable or disable the flat patch trick, use the flags `-flat1` and `-flat2`. For example, to enable it just on the second step, use

    $ nl_bayes noisy.png sigma result.tiff -flat1 0 -flat2 1

The default is `-flat1 1 -flat2 0`.

About this file
---------------

Copyright 2013 IPOL Image Processing On Line http://www.ipol.im/

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
