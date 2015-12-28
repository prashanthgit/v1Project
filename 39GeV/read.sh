#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
      cat "$line"
    done < "$1"
