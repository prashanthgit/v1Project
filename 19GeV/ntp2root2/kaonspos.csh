#!/bin/csh

root.exe -b <<EOF
.L kaonspos.C
kaonspos m("$1")
m.Loop()
.q
EOF

