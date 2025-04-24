# T1MES
Python code to create a T1MES Phantom (See the publication [10.1186/s12968-016-0280-z](https://doi.org/10.1186/s12968-016-0280-z) for more details).
B0 and B1 field are generated with the package [mrtwin](https://infn-mri.github.io/mrtwin/) (see related github repository [here](https://github.com/INFN-MRI/mrtwin))

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
4. Install the project and its dependencies:
   ```bash
   pip install .
   ```

Several new commands will be added to the virtual environment once the installation is completed.
These commands all start with `t1mes_`.

##  ⚙️ Build Static T1MES Phantom
Create a static phantom with the command `t1mes_build_static`.<br>
Here are the available flags: 

| **Flag**                  | **Description** |
|---------------------------|-----------------|
| `-FOV`                    | Field of View of the phantom in millimeters. <br>**Example:** `-FOV 200 200 200` <br>**Default:** `(200, 200, 200)` |
| `-r`, `--resolution`      | Resolution of the phantom (number of voxels per dimension). <br>**Example:** `-r 128 128 5` <br>**Default:** `(128, 128, 5)` |
| `-B0`, `--field_strength` | Static field strength in Tesla. <br>**Example:** `-B0 1.5` <br>**Default:** `1.5` |
| `-PD`                     | Proton density values for `'tubes'` and `'OuterMatrix'`, comma-separated. <br>**Example:** `-PD 1.0,0.5` |
| `-chi`                    | Susceptibility values for `'tubes'` and `'OuterMatrix'`, comma-separated. <br>Use `=` and quotes to avoid issues with negative signs. <br>**Example:** `-chi="-9.,-8."` |
| `-D`                      | Diffusion coefficients for `'tubes'` and `'OuterMatrix'`, comma-separated. <br>**Example:** `-D 0.,0.` |
| `-p`, `--phantom_file`    | Output file path for the generated phantom in `.npz` format. <br>**Example:** `-p phantom.npz` |
| `--plot` / `--no-plot`    | Whether to visualize the generated phantom. <br>Use `--no-plot` to disable plotting. <br>**Default:** `True` |
| `-T2dash`                 | T2′ value (in seconds). <br>**Example:** `-T2dash 0.03` <br>**Default:** `0.03` |
| `--b0_kwargs`             | Additional key-value arguments passed to `mrtwin.b0field()`. <br>**Format:** `key=value` <br>**Example:** `--b0_kwargs offset=10` |
| `--b1_kwargs`             | Additional key-value arguments passed to `mrtwin.b1field()`. <br>**Format:** `key=value` <br>**Example:** `--b1_kwargs scale=0.9` |