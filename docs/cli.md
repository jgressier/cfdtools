# Command Line Tools

For every tool, use option `-h` to get help.

```bash
<command> -h
```

## cfdinfo

`cfdinfo` fully reads all supported formats, converts to an internal mesh and data structure and prints a sum up of available information.

## ic3brief

`ic3brief` currently reads v2 and v3 IC3 files and prints all data headers (mesh and variables) without reading data itself.

## cfdwrite_ic3v2 and cfdwrite_ic3v3

`cfdwrite_ic3v2` and `cfdwrite_ic3v3` specifically write an `IC3 v[23]` file from any supported input.
Some specific options are available to transform mesh or variables:

- `--remove-cell-data varname1 [varname2 ...]` removes the listed names (should be at the end of command line if several names, or followed by `--`)
- `--remove-face-data varname1 [varname2 ...]` removes the listed names (should be at the end of command line if several names, or followed by `--`)
- `--remove-node-data varname1 [varname2 ...]` removes the listed names (should be at the end of command line if several names, or followed by `--`)

`cfdwrite_ic3` is a shortname for last current IC3 writer, namely `cfdwrite_ic3v3`.

## Automatic Format Detection

For all `cfd*` tools, generic file input is supported by an automatic detection of file format through its file extension.
If the extension is missing or unknown, format can be forced with `--fmt <format>`.

- IC3 format with either `.ic3` extension or `--fmt IC3` option. v2 or v3 detection is automatic.
- GMSH format with `.msh` extension or `--fmt GMSH` option. v2.x or v4.x detection is automatic.

## cfdwrite_vtk

## cfdwritecube

## ic3probe_plotline
