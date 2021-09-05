# neuropixel-utils Dev Kit

This directory contains tools for doing development on neuropixel-utils itself. It is for use by the developers of neuropixel-utils, and is not for use by users of it, and does not ship with the project distribution.

neuropixel-utils developers should add this directory to their Matlab path when working on the project.

The whole thing about this directory is that it doesn't care about name collisions in its class and function names; it assumes that neuropixel-utils is the only thing you're working on in your session. It names functions whatever it wants, and doesn't put them inside packages. (This is intentional, to make it convenient for neuropixel-utils developers to call them without having to type out the package name all the time.) So it is _not_ suitable for redistribution to users.
