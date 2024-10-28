# Command Line Tools

For every tool, use option `-h` to get help.

```bash
<command> -h
```


## generic writers

Some specific options are available to transform mesh or variables when writing a new mesh/solution file:

- `--remove-cell-data varname1 [varname2 ...]` removes the listed names (should be at the end of command line if several names, or followed by `--`)
- `--remove-node-data varname1 [varname2 ...]` removes the listed names (should be at the end of command line if several names, or followed by `--`)


::: cfdtools._cli
    options:
        show_signature: false
        docstring_section_style: list
        show_docstring_parameters: false
        show_docstring_returns: false
        show_source: false
        heading_level: 2
        parameter_headings: false
        show_root_toc_entry: false
        filters:
            - "!.*cli.*"
