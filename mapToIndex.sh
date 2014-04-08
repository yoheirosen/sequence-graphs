#!/usr/bin/env bash
# Script to run the "mapToIndex" project

# Get the current directory (see <http://stackoverflow.com/a/246128/402891>)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Run the command. Make sure to build with "sbt stage" first.
$DIR/mapToIndex/target/universal/stage/bin/maptoindex "$@"
