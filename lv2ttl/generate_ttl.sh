#!/bin/bash
LIB_EXT=$1
rm manifest.ttl || true
sed "s/@LIB_EXT@/$LIB_EXT/" lv2ttl/manifest.ttl.in > lv2ttl/manifest.ttl
