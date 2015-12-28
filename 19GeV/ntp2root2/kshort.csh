#!/bin/csh

root.exe -b <<EOF
.L kshort.C
kshort m("$1")
m.Loop()
.q
EOF

