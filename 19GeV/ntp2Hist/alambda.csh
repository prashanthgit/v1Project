#!/bin/csh

root.exe -b <<EOF
.L alambda.C
alambda m
m.Loop()
.q
EOF

