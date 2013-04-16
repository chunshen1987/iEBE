#! /usr/bin/env bash

./iS.e

./resonance.e

./iInteSp.e

(cd extractThermal; ./extractThermal3.sh ../results)
