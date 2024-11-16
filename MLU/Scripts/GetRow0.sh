#!/usr/bin/env bash

awk -f ${0%.*}.awk -v field="${field:-corr}" "$1"
