# protein_analysis
Protein analysis project.

Instructions
------------

1. Clone repository:
```
$ git clone https://github.com/JuanjoAlegria/protein_analysis.git
```

2. Create the environment:
```
$ conda create --name protein_analysis --file env-explicit.txt
```

3. Activate the environment
```
$ conda activate protein_analysis
```

4. Execute the script

```
$ python -m src.protein_analysis --df_features_path path/to/features.csv --save_dir path/to/folder
```

Other parameters
----------------

- --proteins \<proteins\>: List of proteins that the script should consider. All the other proteins will not be taken into account. Example:

```
$ python -m src.protein_analysis
--df_features_path path/to/features.csv
--save_dir path/to/folder
--proteins ABI2_A3_NO PRKACA_F2_NO SHC1_G1_NO
```

- --no_labels: Flag that indicates if the plot should include the name of each protein. Example
```
$ python -m src.protein_analysis
--df_features_path path/to/features.csv
--save_dir path/to/folder
--no_labels
```

- --min_x \<min_x\>: Minimum X coordinate in the plot. If not included, it will be infered from the samples.
- --max_x \<max_x\>: Maximum X coordinate in the plot. If not included, it will be infered from the samples.
- --min_y \<min_y\>: Minimum Y coordinate in the plot. If not included, it will be infered from the samples.
- --max_y \<max_y\>: Maximum Y coordinate in the plot. If not included, it will be infered from the samples.



