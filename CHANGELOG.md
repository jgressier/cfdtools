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

### [0.x.x](https://pypi.org/project/cfdtools/) (2021-xx-xx)

#### changed

- avoid warnings with `vanalbada` and `vanleer` limiters when uniform flows
- analytical 1D solution for nozzle flows in `solution.euler_nozzle`
- improve test coverage
- optimize some mesh computation

#### fixed

