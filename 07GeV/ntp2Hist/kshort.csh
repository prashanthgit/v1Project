#!/bin/csh

root.exe -b <<EOF
.L kshort.C
kshort m
m.Loop()
.q
EOF

