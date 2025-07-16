#!/bin/bash

clang main.c cocoa_bridge.m -framework Cocoa -o cube

./cube