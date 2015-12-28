#!/bin/csh

root.exe -b <<EOF
.L aproton.C
aproton m("$1")
m.Loop()
.q
EOF

