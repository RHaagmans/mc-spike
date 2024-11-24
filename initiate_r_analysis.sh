#!/bin/bash

help_message()
{
   # Display Help
   echo "Process data from the MC Spike nextflow pipeline and then generate reports."
   echo "Optionally provide the path to a <data_folder> containing the output from the"
   echo "pipeline"
   echo
   echo "Usage:"
   echo "initiate_r_analysis.sh [data_folder] [-h]"
   echo 
   echo "Options:"
   echo "-h     Print this help message."
   echo
}

while getopts ":h" option; do
   case $option in
      h)
         help_message
         exit;;
     \?) # Invalid option
         echo "Error: Invalid option"
         help_message
         exit;;
   esac
done

echo
echo "Processing data..."
echo

Rscript "R/process_data.R" $1

echo
echo "Generating reports..."
echo

Rscript "generate_reports.R" $1
