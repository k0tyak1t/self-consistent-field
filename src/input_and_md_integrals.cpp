#include "input_and_md_integrals.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "matrix.h"
#include "mur_dav_prim_integrals.h"
#include "standard_matrices.h"
#include "utils.h"
#define ANGSTROM_TO_BOHR 1.8897259886

bool InputAndMDIntegrals::readGeom(char *filename) {
  std::ifstream inp(filename);
  int num_atoms;
  inp >> num_atoms;
  std::string cur_str;
  getline(inp, cur_str); // skipping of comment line
  getline(inp, cur_str);
  Atom cur_atom;
  while (inp >> cur_str >> cur_atom.x >> cur_atom.y >> cur_atom.z) {
    cur_atom.x *= ANGSTROM_TO_BOHR;
    cur_atom.y *= ANGSTROM_TO_BOHR;
    cur_atom.z *= ANGSTROM_TO_BOHR;
    cur_atom.q = get_nuclear_charge(cur_str);
    if (cur_atom.q == 0) {
      throw std::runtime_error("Incorrect atom.");
    }
    atoms.push_back(cur_atom);
  }
  inp.close();
  if (num_atoms != int(atoms.size())) {
    throw std::runtime_error("Incorrect number of atoms.");
  }
  return true;
}

bool InputAndMDIntegrals::readBasisLib(char *filename) {
  std::ifstream inp(filename);
  std::string tmp_str;
  double C, a;
  int nFunc, tmp_int, idAtom;
  vector<pair<int, vector<pair<double, double>>>> singleAtomBL;
  vector<pair<double, double>> singleBF;

  inp >> tmp_str;

  while (true) {
    idAtom = get_nuclear_charge(tmp_str);
    if (idAtom < 1)
      return false;
    singleAtomBL.clear();
    inp >> tmp_str;
    while (get_orbital_momentum(tmp_str) >= 0) {
      inp >> nFunc;
      singleBF.clear();
      for (int i = 0; i < nFunc; i++) {
        inp >> tmp_int >> a >> C;
        singleBF.push_back(make_pair(a, C));
      }
      singleAtomBL.push_back(
          make_pair(get_orbital_momentum(tmp_str), singleBF));
      if (!(inp >> tmp_str)) {
        inp.close();
        basisLib[idAtom] = singleAtomBL;
        return true;
      }
    }
    basisLib[idAtom] = singleAtomBL;
  }
}

bool InputAndMDIntegrals::calc(StandardMatrices &std_m) {
  std_m.init(basisFunctions.size());
  //	Вычисление энергии электрон-электронного взаимодействия и числа
  // электронов из условия электронейтральности
  double total_Vnn = 0;
  int num_el = 0;
  for (int i = 0; i < int(atoms.size()); i++) {
    num_el += atoms[i].q;
    for (int j = i + 1; j < int(atoms.size()); j++)
      total_Vnn += atoms[i].q * atoms[j].q /
                   sqrt(pow(atoms[i].x - atoms[j].x, 2) +
                        pow(atoms[i].y - atoms[j].y, 2) +
                        pow(atoms[i].z - atoms[j].z, 2));
  }
  std_m.set_total_Vnn(total_Vnn);
  std_m.set_num_el(num_el);

  //	расчет матричных элементов S,T,Hcore
  int ii, jj, kk, ll;
  ii = jj = 0;
  for (vector<vector<pair<int, SingleBasisFunction>>>::iterator i =
           basisFunctions.begin();
       i != basisFunctions.end(); i++) {
    jj = 0;
    for (vector<vector<pair<int, SingleBasisFunction>>>::iterator j =
             basisFunctions.begin();
         j != basisFunctions.end(); j++) {
      double Sij = 0, Tij = 0, Vij = 0;
      for (vector<pair<int, SingleBasisFunction>>::iterator i1 = i->begin();
           i1 != i->end(); i1++)
        for (vector<pair<int, SingleBasisFunction>>::iterator j1 = j->begin();
             j1 != j->end(); j1++) {
          double Sij_t = 0, Tij_t = 0, Vij_t = 0;
          SingleBasisFunction f1 = i1->second;
          SingleBasisFunction f2 = j1->second;
          for (vector<pair<double, double>>::iterator i2 = f1.ai_ci.begin();
               i2 != f1.ai_ci.end(); i2++)
            for (vector<pair<double, double>>::iterator j2 = f2.ai_ci.begin();
                 j2 != f2.ai_ci.end(); j2++) {
              double C = i2->second * j2->second;
              Sij_t += MurDavPrimIntegrals().Sij(
                           f1.nx, f1.ny, f1.nz, f1.x, f1.y, f1.z, i2->first,
                           f2.nx, f2.ny, f2.nz, f2.x, f2.y, f2.z, j2->first) *
                       C;
              Tij_t += MurDavPrimIntegrals().Tij(
                           f1.nx, f1.ny, f1.nz, f1.x, f1.y, f1.z, i2->first,
                           f2.nx, f2.ny, f2.nz, f2.x, f2.y, f2.z, j2->first) *
                       C;
              for (int k = 0; k < int(atoms.size()); k++)
                Vij_t += MurDavPrimIntegrals().Vij(
                             f1.nx, f1.ny, f1.nz, f1.x, f1.y, f1.z, i2->first,
                             f2.nx, f2.ny, f2.nz, f2.x, f2.y, f2.z, j2->first,
                             atoms[k].x, atoms[k].y, atoms[k].z) *
                         C * atoms[k].q;
            }
          Sij += Sij_t * i1->first * j1->first;
          Tij += Tij_t * i1->first * j1->first;
          Vij += Vij_t * i1->first * j1->first;
        }

      std_m.S[ii][jj] = Sij;
      std_m.T[ii][jj] = Tij;
      std_m.Ven[ii][jj] = -Vij;
      std_m.H[ii][jj] = Tij - Vij;
      jj++;
    }
    ii++;
  }
  //	расчет двуэлектронных интегралов
  ii = jj = kk = ll = 0;
  const int ii_end = int(basisFunctions.size());

  unsigned char progress;

  for (vector<vector<pair<int, SingleBasisFunction>>>::iterator i =
           basisFunctions.begin();
       i != basisFunctions.end(); i++) {
    jj = kk = ll = 0;
    for (vector<vector<pair<int, SingleBasisFunction>>>::iterator j =
             basisFunctions.begin();
         j != basisFunctions.end(); j++) {
      kk = ll = 0;
      for (vector<vector<pair<int, SingleBasisFunction>>>::iterator k =
               basisFunctions.begin();
           k != basisFunctions.end(); k++) {
        ll = 0;
        for (vector<vector<pair<int, SingleBasisFunction>>>::iterator l =
                 basisFunctions.begin();
             l != basisFunctions.end(); l++) {
          double Vijkl = 0;
          for (vector<pair<int, SingleBasisFunction>>::iterator i1 = i->begin();
               i1 != i->end(); i1++)
            for (vector<pair<int, SingleBasisFunction>>::iterator j1 =
                     j->begin();
                 j1 != j->end(); j1++)
              for (vector<pair<int, SingleBasisFunction>>::iterator k1 =
                       k->begin();
                   k1 != k->end(); k1++)
                for (vector<pair<int, SingleBasisFunction>>::iterator l1 =
                         (*l).begin();
                     l1 != (*l).end(); l1++) {
                  double Vijkl_t = 0;
                  SingleBasisFunction f1 = i1->second;
                  SingleBasisFunction f2 = j1->second;
                  SingleBasisFunction f3 = k1->second;
                  SingleBasisFunction f4 = l1->second;
                  for (vector<pair<double, double>>::iterator i2 =
                           f1.ai_ci.begin();
                       i2 != f1.ai_ci.end(); i2++)
                    for (vector<pair<double, double>>::iterator j2 =
                             f2.ai_ci.begin();
                         j2 != f2.ai_ci.end(); j2++)
                      for (vector<pair<double, double>>::iterator k2 =
                               f3.ai_ci.begin();
                           k2 != f3.ai_ci.end(); k2++)
                        for (vector<pair<double, double>>::iterator l2 =
                                 f4.ai_ci.begin();
                             l2 != f4.ai_ci.end(); l2++)
                          Vijkl_t +=
                              MurDavPrimIntegrals().Vijkl(
                                  f1.nx, f1.ny, f1.nz, f1.x, f1.y, f1.z,
                                  i2->first, f2.nx, f2.ny, f2.nz, f2.x, f2.y,
                                  f2.z, j2->first, f3.nx, f3.ny, f3.nz, f3.x,
                                  f3.y, f3.z, k2->first, f4.nx, f4.ny, f4.nz,
                                  f4.x, f4.y, f4.z, l2->first) *
                              i2->second * j2->second * k2->second * l2->second;
                  Vijkl +=
                      Vijkl_t * i1->first * j1->first * k1->first * l1->first;
                }
          std_m.set_eri(ii, jj, kk, ll, Vijkl);
          ll++;
        }
        kk++;
      }
      jj++;
    }
    ii++;
    progress = 100 * ii / ii_end;
    display_progress(progress, "ERI calculation: ");
  }
  std::cout << std::endl;
  //	Матрица X
  int nAO = basisFunctions.size();
  std_m.X.init(nAO);
  double *evec = new double[nAO * nAO];
  double *eval = new double[nAO];
  std_m.S.eigen_vv(evec, eval);
  matrix s;
  s.init(nAO);
  for (int i = 0; i < nAO; i++) { // создание матрицы S^(-0.5)
    for (int j = 0; j < nAO; j++) {
      s[i][j] = 0;
    }
    s[i][i] = 1.0 / sqrt(eval[i]);
  }
  delete[] eval;
  matrix U;
  U.init(nAO);
  U.from_array(evec);
  delete[] evec;
  matrix vs;
  vs.init(nAO);
  std_m.X = U.T() * s * U;
  return true;
}

void InputAndMDIntegrals::printBasis() {
  int nBFs = 0;
  for (vector<vector<pair<int, SingleBasisFunction>>>::iterator it =
           basisFunctions.begin();
       it != basisFunctions.end(); it++) {
    std::cout << "Basis function #" << ++nBFs << ":\n";
    for (vector<pair<int, SingleBasisFunction>>::iterator jt = it->begin();
         jt != it->end(); jt++) {
      std::cout << "  " << jt->first;
      std::cout << "  X=" << jt->second.x;
      std::cout << "  Y=" << jt->second.y;
      std::cout << "  Z=" << jt->second.z;
      std::cout << "  nX=" << jt->second.nx;
      std::cout << "  nY=" << jt->second.ny;
      std::cout << "  nZ=" << jt->second.nz << '\n';
      for (vector<pair<double, double>>::iterator kt = jt->second.ai_ci.begin();
           kt != jt->second.ai_ci.end(); kt++)
        std::cout << "         " << (*kt).first << "  " << (*kt).second << '\n';
    }
  }
}

bool InputAndMDIntegrals::set_c2v_z(int *symmAO) {
  int nBFs = 0;
  std::cout << "Проверка и установка C2v симметрии для линейных молекул, "
               "ориентированных на оси Z\n";
  bool isGood = true;
  for (vector<vector<pair<int, SingleBasisFunction>>>::iterator it =
           basisFunctions.begin();
       it != basisFunctions.end(); it++) {
    std::cout << "Basis function #" << ++nBFs << ":\n";
    int symm = 0;
    for (vector<pair<int, SingleBasisFunction>>::iterator jt = it->begin();
         jt != it->end(); jt++) {
      std::cout << "  " << jt->first;
      std::cout << "  X=" << jt->second.x;
      std::cout << "  Y=" << jt->second.y;
      std::cout << "  Z=" << jt->second.z;
      std::cout << "  nX=" << jt->second.nx;
      std::cout << "  nY=" << jt->second.ny;
      std::cout << "  nZ=" << jt->second.nz << '\n';
      if ((jt->second.x != 0) || ((*jt).second.y != 0)) {
        std::cout << "Эта функция не лежит на оси z, симметрия не будет "
                     "установлена\n";
        isGood = false;
      }
      int symm_ = 2 * (jt->second.ny % 2) + ((*jt).second.nx % 2) + 1;
      if (symm == 0)
        symm = symm_;
      if (symm != symm_)
        symm = -1;
    }
    if ((!isGood) || (symm < 1)) {
      std::cout
          << "Эта функция несимметрична, симметрия не будет установлена\n";
      return false;
    }
    std::cout << "Function symmetry: ";
    if (symm == 1) {
      std::cout << "A1\n";
    }

    if (symm == 2) {
      std::cout << "B1\n";
    }

    if (symm == 3) {
      std::cout << "B2\n";
    }

    if (symm == 4) {
      std::cout << "A2\n";
    }

    symmAO[nBFs - 1] = symm;
  }

  return true;
}

bool InputAndMDIntegrals::makeBasis() {
  int i = 0;
  int progress = 0;
  std::string leading_str = "Loading basis: ";
  vector<pair<double, double>> kk1;
  vector<pair<double, double>> kk2;
  kk1 = kk2;
  std::cout << "\n";
  for (vector<Atom>::iterator atom = atoms.begin(); atom != atoms.end();
       ++atom) {
    display_progress(progress, leading_str);
    i++;
    progress = 100 * (i + 1) / atoms.size();

    map<int, vector<pair<int, vector<pair<double, double>>>>>::iterator basis =
        basisLib.find(atom->q);
    if (basis == basisLib.end()) {
      throw std::runtime_error("Failed to load atom with charge: " +
                               std::to_string(atom->q));
    }

    for (vector<pair<int, vector<pair<double, double>>>>::iterator it1 =
             basis->second.begin();
         it1 != basis->second.end(); ++it1) {
      int l = it1->first;
      SingleBasisFunction sbf_;
      vector<pair<int, SingleBasisFunction>> mbf_;
      sbf_.ai_ci = it1->second;
      if (l != 2)
        for (int lx = 0; lx <= l; lx++)
          for (int ly = 0; ly <= l - lx; ly++) {
            mbf_.clear();
            sbf_.setLC(lx, ly, l - lx - ly, atom->x, atom->y, atom->z);
            mbf_.push_back(make_pair(1, sbf_));
            basisFunctions.push_back(mbf_);
          }
      // cлучай чистых d
      if (l == 2) {
        //      dxy
        mbf_.clear();
        sbf_.setLC(1, 1, 0, atom->x, atom->y, atom->z);
        sbf_.ai_ci = it1->second;
        mbf_.push_back(make_pair(1, sbf_));
        basisFunctions.push_back(mbf_);
        //      dxz
        mbf_.clear();
        sbf_.setLC(1, 0, 1, atom->x, atom->y, atom->z);
        mbf_.push_back(make_pair(1, sbf_));
        basisFunctions.push_back(mbf_);
        //      dyz
        mbf_.clear();
        sbf_.setLC(0, 1, 1, atom->x, atom->y, atom->z);
        mbf_.push_back(make_pair(1, sbf_));
        basisFunctions.push_back(mbf_);

        //      d(x2-y2)
        mbf_.clear();
        sbf_.setLC(2, 0, 0, atom->x, atom->y, atom->z);
        mbf_.push_back(make_pair(1, sbf_));
        sbf_.setLC(0, 2, 0, atom->x, atom->y, atom->z);
        mbf_.push_back(make_pair(-1, sbf_));
        basisFunctions.push_back(mbf_);
        //	d(2z2-x2-y2)
        mbf_.clear();
        sbf_.setLC(0, 0, 2, atom->x, atom->y, atom->z);
        mbf_.push_back(make_pair(2, sbf_));
        sbf_.setLC(2, 0, 0, atom->x, atom->y, atom->z);
        mbf_.push_back(make_pair(-1, sbf_));
        sbf_.setLC(0, 2, 0, atom->x, atom->y, atom->z);
        mbf_.push_back(make_pair(-1, sbf_));
        basisFunctions.push_back(mbf_);
      }
    }
  }

  for (vector<vector<pair<int, SingleBasisFunction>>>::iterator it =
           basisFunctions.begin();
       it != basisFunctions.end(); it++)
    for (vector<pair<int, SingleBasisFunction>>::iterator jt = it->begin();
         jt != it->end(); jt++)
      jt->second.renorm_ci();
  return true;
}

int InputAndMDIntegrals::get_nuclear_charge(const std::string &str_) const {
  std::string str(str_);
  std::transform(str.begin(), str.end(), str.begin(), toupper);
  if ((str == "H") || (str == "HYDROGEN"))
    return 1;
  if ((str == "HE") || (str == "HELIUM"))
    return 2;
  if ((str == "LI") || (str == "LITHIUM"))
    return 3;
  if ((str == "BE") || (str == "BERYLLIUM"))
    return 4;
  if ((str == "B") || (str == "BORON"))
    return 5;
  if ((str == "C") || (str == "CARBON"))
    return 6;
  if ((str == "N") || (str == "NITROGEN"))
    return 7;
  if ((str == "O") || (str == "OXYGEN"))
    return 8;
  if ((str == "F") || (str == "FLUORINE"))
    return 9;
  if ((str == "AR") || (str == "ARGON"))
    return 18;

  throw std::runtime_error("Atom were not recognized! Incorrect input");
  return 0;
}

int InputAndMDIntegrals::get_orbital_momentum(const std::string &str_) const {
  std::string str(str_);
  std::transform(str.begin(), str.end(), str.begin(), toupper);
  if ((str == "S") || (str == "L0"))
    return 0;
  if ((str == "P") || (str == "L1"))
    return 1;
  if ((str == "D") || (str == "L2"))
    return 2;
  if ((str == "F") || (str == "L3"))
    return 3;
  if ((str == "G") || (str == "L4"))
    return 4;
  if ((str == "H") || (str == "L5"))
    return 5;
  if ((str == "I") || (str == "L6"))
    return 6;
  if ((str == "J") || (str == "L7"))
    return 7;
  if ((str == "K") || (str == "L8"))
    return 8;
  if (str == "L9")
    return 9;
  if (str == "L10")
    return 10;
  if (str == "L11")
    return 11;
  return -1;
}
