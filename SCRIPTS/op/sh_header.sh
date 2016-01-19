#!/bin/bash
vim_opts="-n -p"
this_script=` basename $0 `
me=$this_script
this_script_dir="$shd"
root_dir="$shd/../"
data_dir="$root_dir/data/"

base_output_dir="out"
output_dir="$data_dir/$base_output_dir"

source $shd/f.sh

