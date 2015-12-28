#!/bin/csh

root.exe -b <<EOF
.L proton.C
proton m("$1")
m.Loop()
.q
EOF

