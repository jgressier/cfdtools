# IC3 format V2 and V3

## restart header

- `UGP_IO_MAGIC_NUMBER` to test byte swap (little endian, big endian)
- version

## section headers

order of sections is arbitrary and can be one of the following

- `UGP_IO_NO_FA_CV_NOOFA_COUNTS`: first integers are
  1. no_count
  2. fa_count
  3. cv_count
  4. noofa_count
- `UGP_IO_NO_CHECK`
- `UGP_IO_FA_CHECK`
- `UGP_IO_CV_CHECK`
- `UGP_IO_NOOFA_I_AND_V`: CSR connectivity of face2node
- `UGP_IO_CVOFA`
- `UGP_IO_FA_ZONE`
- `UGP_IO_CV_PART`
- `UGP_IO_X_NO`
- `UGP_IO_DATA`
- `UGP_IO_I0`
- `UGP_IO_D0`
- `UGP_IO_FA_II1`
- `UGP_IO_FA_D1`
- `UGP_IO_FA_D3`
- `UGP_IO_CV_II1`
- `UGP_IO_CV_D1`
- `UGP_IO_CV_D3`
- `UGP_IO_CV_D33`: cell based and double type 3x3 tensor data
- `UGP_IO_NO_II1`: node based and int scalar data
- `UGP_IO_NO_D1`: node based and double scalar data
- `UGP_IO_NO_D3`: node based and double vector data
- `UGP_IO_EOF`: end of file

## mesh data

order of sections is arbitrary

### `NOOFA` face2node connectivity

- header with number of nodes, faces and cells

## parameters
## variables and fieds

Though the format can be read in any order. The current (v2 and v3) order in IC3 code is

- scalars
  - double float: nodes, faces, and cells
  - long int: nodes, faces, and cells
- vectors
  - double float: nodes, faces, and cells
- tensors
  - double float: cells only
