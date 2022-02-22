#!/bin/bash

g++ -O2 -pthread -march=native -funroll-loops planar_monomials.c -o planar_monomials -lntl -lgmp -lm


