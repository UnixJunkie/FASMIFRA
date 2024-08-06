#!/usr/bin/env bash

sed -E 's/\[\*:[0-9]+\]\[\*:[0-9]+\]//g' $1
