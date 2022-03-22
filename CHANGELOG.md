# Changelog & Release Notes

## Upgrading

To upgrade to the latest version of `cfdtools` use `pip`:

```bash
pip install cfdtools --upgrade
```

You can determine your currently installed version using this command:

```bash
pip show cfdtools
```

## Versions

### [0.2.0](https://pypi.org/project/cfdtools/) (2022-03-22)

#### new

- command line `ic3brief` only section headers
- command line writers can remove date with `--remove-[cell/node]-data` options

#### changed

- ensure safe new names for writers to avoid deletion of read files

### [0.1.0](https://pypi.org/project/cfdtools/) (2022-03-11)

#### new

- command line `cfdwrite_ic3v2`, `cfdwrite_ic3v3`
- new `v3 IC3` writer

#### changed

- improve IC3 v2 and v3 with parameters and actual partition

### [0.0.2](https://pypi.org/project/cfdtools/) (2022-03-04)

#### new

- command line `cfdinfo`
- new `v2 IC3` writer
#### fixed

- v2 IC3 reader 
