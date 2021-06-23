
import psi4
psi4.geometry("""
Ne
Ne 1 3.0
""")
psi4.set_options({"freeze_core": "True"})
psi4.energy("ccsd(t)/cc-pvtz")
