#!/bin/bash

autoreconf --install
automake --add-missing --copy &> /dev/null
exit 0

