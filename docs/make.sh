#!/bin/bash
# Idea of this script taken from: http://steven.casagrande.io/articles/travis-ci-and-if-statements/

set -ev

# Build the doc
if [ "$TRAVIS_OS_NAME" == "linux" ]; then
    julia -e 'Pkg.add("PhyloPlots")';
    julia -e 'using Pkg; ps=PackageSpec(name="Documenter", version="0.19"); Pkg.add(ps); Pkg.pin(ps)';
    julia -e 'cd(Pkg.dir("PhyloNetworks")); include(joinpath("docs", "make.jl"))';
fi

exit 0;
