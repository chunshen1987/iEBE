#! /usr/bin/env bash

# resonance calculation
cp resonance8/reso.inp ./
./resonance.e
rm reso.inp

# calculate flow
./iInteSp.e

