#!/bin/sh

VERSION=${1:-latest}
docker build -t jkadbear/gr-lora:$VERSION .
