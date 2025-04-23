# T1MES
Python code to create a T1MES Phantom (See the publication [10.1186/s12968-016-0280-z](https://doi.org/10.1186/s12968-016-0280-z) for more details).

## Install

1. Download the repository:
   ```bash
   # clone project
   git clone https://github.com/CyprienBouton/T1MES
   cd T1MES
   ```
2. Create a virtual environment (Conda is strongly recommended):
   ```bash
   # create conda environment
   conda create -n T1MES python=3.10
   conda activate T1MES
   ```
4. Install the project in editable mode and its dependencies:
   ```bash
   pip install .
   ```

Several new commands will be added to the virtual environment once the installation is completed.
These commands all start with `t1mes_`.