# GEOP 603 term project: Fourier finite-difference wave propagation

Replication and extension of Song, X., and S. Fomel (2011), Fourier finite-difference wave propagation, *Geophysics* 76(5), T123–T129.

- `notebooks/`: implementation, digitization and verification notebooks
- `report/`: term report and talk material

## Reproducing (C1 verification)

    pip install numpy scipy matplotlib jupyter pytest
    pytest -q notebooks/tests
    cd notebooks && jupyter nbconvert --to notebook --execute --inplace C1_verification.ipynb

Figures are written to `report/figures/`.
