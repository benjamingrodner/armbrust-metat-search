#!/bin/bash
# Input $0 is the path to your etetoolkit directory
if [ -d "$HOME" ]; then
  cd $HOME
  ln -sf -T $0 .etetoolkit
else 
  printf '%s\n' "Can't access home directory. Must mount home directory or specify a scratch directory (e.g. apptainer run --no-home -S $HOME)" >&2  # write error message to stderr
  exit 1                                  # or exit $test_result
fi
