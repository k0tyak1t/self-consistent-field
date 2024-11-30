#include "output.h"

output::output() { }
output::~output() { }

int output::run(const MOs& mos) const
{
    std::cout << "\n--- Orbitals ---" << std::endl;
    mos.C.print();
    std::cout << "\n--- Final energies ---" << std::endl;
    int nAO = mos.C.get_size();
    for (int i = 0; i < nAO; i++)
        std::cout << std::setw(8) << std::setprecision(8) << mos.get_mo_energy(i) << ' ';
    std::cout << '\n';
    std::cout << "\n total_energy=" << mos.get_total_energy() << '\n';
    return 0;
}

void output::print_energies_only(const MOs& mos) const
{
    std::cout << "\n--- Final energies ---" << std::endl;
    int nAO = mos.C.get_size();
    for (int i = 0; i < nAO; i++) {
        std::cout << '#' << i + 1 << ':' << std::setw(8) << std::setprecision(8) << mos.get_mo_energy(i) << std::endl;
    }
    std::cout << "Total energy: " << mos.get_total_energy() << std::endl;
}
