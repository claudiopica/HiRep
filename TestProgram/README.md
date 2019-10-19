## HiRep tests

This directory contains test programs for the HiRep code.

To run all the available tests using the current `MkFlags` settings, execute

```
make tests && make runalltests
```
from the root TestProgram folder.

Test are organized in modules, each with its own folder. To run test for a single module with the current `MkFlags` settings, run:

```
make && make runtests
```
from a module folder.

All tests within a module must be able to run without any command line arguments and must return `0` on success. Any input file necessary for the run must be available to the test. The default input file is called `input_file`.

For tests compiled with mpi support, a simple wrapper is provided to run the test, called `mpi_test_wrapper.pl`. The wrapper reads the default `input_file` to determine the number of mpi processes requested and run the test via `mpirun`. 


## CI Pipeline

Tests are executed automatically on each commit by using Github Actions. Two workflows are defined: one for non-mpi (serial) compilation and one for mpi (parallel). The two workflows are otherwise identical.
Each CI workflow define a matrix on compilation flags (number of colors, representation, etc).

The script _run_tests.sh_ automate the process of: creating the appropriate `MkFlags`; compiling the tests; and running all the test. Run:
```
run_tests.sh --help
```
for a list of available options and their default values. For example:
```
run_tests.sh -mpi -n 2 -r FUND
```
will write a new `MkFlags` file for 2 colors and fermions in the fundamental representation (for all other required flags default values will be used).

Tests on Github Actions are executed inside a Docker container. The Dockerfile for the container is in the TestProgram directory, and the image is hosted at `docker://cpica/hr-tests`.

To reproduce the execution in the docker container locally, run:

```
export GITHUB_WORKSPACE=/github/workspace
docker run --workdir /github/workspace -e GITHUB_WORKSPACE -v $(pwd)/HiRep:/github/workspace cpica/hr-tests ./TestProgram/run_tests.sh -mpi -n 2 -r FUND
```
e.g. for the case of 2 gauge colors and fermions in the fundamental representation.

Github Actions workflows are defined in the `.github/` folder at the root of the `HiRep` repository. 





