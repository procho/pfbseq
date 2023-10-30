# Setting up conda enviroment for scanpy and jupyter

## Create project environment

```bash
conda create -n final_project python=3.8
```

## Activate created environment

```bash
conda activate final_project
```

## Install packages

```bash
conda install -c conda-forge jupyter numpy pandas scanpy matplotlib
```

## May need to update matplotlib

```bash
conda update matplotlib
``` 

## Open jupyter notebook 

```bash
jupyter notebook
```
