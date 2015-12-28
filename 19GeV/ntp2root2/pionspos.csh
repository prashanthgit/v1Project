#!/bin/csh

root.exe -b <<EOF
.L pionspos.C
pionspos m("$1")
m.Loop()
.q
EOF

