#!/bin/bash

clang -framework Cocoa -framework Metal -framework MetalKit -framework QuartzCore -framework Foundation main.c metal_bridge.m -o cube

./cube