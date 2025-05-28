#include "preset.h"

int main()
{
    input Input{"INPUT"};
    cpi_rand &Rand = cpi_rand::get_instance();
    box Box{Input};
    output Output{Input, Box};
    simu Simulation{Input, Output, Rand, Box};
    Input.print("ALL_INPUT");
    Simulation.evolve();
    return 0;
}