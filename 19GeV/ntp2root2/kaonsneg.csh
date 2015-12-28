#!/bin/csh

root.exe -b <<EOF
.L kaonsneg.C
kaonsneg m("$1")
m.Loop()
.q
EOF

