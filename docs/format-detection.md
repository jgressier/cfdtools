# Automatic Format Detection

For all `cfd*` tools, generic file input is supported by an automatic detection of file format through its file extension.
If the extension is missing or unknown, format can be forced with `--fmt <format>`.

- IC3 format with either `.ic3` extension or `--fmt IC3` option. v2 or v3 detection is automatic.
- GMSH format with `.msh` extension or `--fmt GMSH` option. v2.x or v4.x detection is automatic.