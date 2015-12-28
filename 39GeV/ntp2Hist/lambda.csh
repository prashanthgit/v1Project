#!/bin/csh

root.exe -b <<EOF
.L lambda.C
lambda m
m.Loop()
.q
EOF

