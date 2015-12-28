#!/bin/csh

root.exe -b <<EOF
.L alambda.C
alambda m("$1")
m.Loop()
.q
EOF

