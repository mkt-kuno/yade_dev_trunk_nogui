#!/bin/bash

set -e
for i in bullseye focal jammy bookworm trixie noble forky resolute
do
    aptly repo create -distribution=$i -component=main yadedaily-$i
done
