import numpy as np

# elemental abundances
ele = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At SiO"
aele = ele.split()
nele = len(aele)


def main():
    rand_abun = np.random.rand(nele)
    norm_abun = rand_abun / np.sum(rand_abun)

    # grid sizes
    nt = 50
    ng = 32

    # grid states
    t0 = 0.0
    tf = 1.0e4
    temp0 = 1.0e4
    tempf = 1.0e2
    vol0 = 1.0e22
    rho0 = 1.0e-2
    prs0 = 1.0e3

    dt = (tf - t0) / float((nt - 1))
    dtemp = (tempf - temp0) / float((nt - 1))

    with open("bm_abundances.dat", "w") as af:
        print(f"grid {ele}", file=af)
        for g in range(ng):
            row = f"{g + 1}"
            for e in range(nele):
                row += f" {norm_abun[e]}"
            print(f"{row}", file=af)

    with open("bm_traj.dat", "w") as tf:
        for i in range(nt):
            tnow = t0 + float(i) * dt
            tempnow = temp0 + float(i) * dtemp
            print(f"{tnow}", file=tf)
            for g in range(ng):
                print(f"{g + 1} {tempnow} {vol0} {rho0} {prs0} 0 0", file=tf)


if __name__ == "__main__":
    main()
