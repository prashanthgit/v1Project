#!/bin/csh

root.exe -b <<EOF
.L pionsneg.C
pionsneg m("$1")
m.Loop()
.q
EOF

