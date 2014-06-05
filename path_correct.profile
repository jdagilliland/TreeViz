#!/bin/sh
if [ -d "/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin" ] ; then
    PATH="/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin:$PATH"
else
    echo "I could not find the MacPorts version of python2.7; is it installed?"
fi
