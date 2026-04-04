# Oscar

Oscar is a set of function to read Finite Element output files in Hdf5 format

## installation

./install

## Examples

The `h5tovtk` command accepts input filenames with or without the `.h5`
extension. For example, both `h5tovtk simulation` and
`h5tovtk simulation.h5` work for single-file results.

## License

[LICENSE.txt](LICENSE.txt)

## Develoment log

1.1 Fully functional

1.2 Added support for reduced files (for constant meshes).

1.2.1 Improved speed by creating elemNodes array first.

