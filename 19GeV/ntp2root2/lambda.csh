#!/bin/csh

root.exe -b <<EOF
.L lambda.C
lambda m("$1")
m.Loop()
.q
EOF

