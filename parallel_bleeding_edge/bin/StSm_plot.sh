#!/bin/sh

if ! command -v splash &> /dev/null; then
    echo "splash is not installed."
    echo "Please refer to https://splash-viz.readthedocs.io/en/latest/getting-started.html for installation instructions."
    echo "Your operating system is: $(uname -s)"
    exit 1
else
    SPLASH=$(which splash)
    $SPLASH -f starsmasher "$@"
fi
