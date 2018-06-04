# Sonography Interpolation

## Building

## Unix-like

```bash
mkdir -p build
cd build
cmake ..
make
```

The executable `usi_interpolate` will be generated under "`build`".

## Windows

```cmd
mkdir build
cd build
cmake ..
msbuild usi_interpolation.sln
```
The executable `usi_interpolate.exe` will be generated under "`build/Debug`".

## Usage

    usi_interpolate FILE

e.g.

    usi_interpolate image.dat

Then `raw.png` and `interpolated.png` will be created.

## Sample outputs

Example outputs are under the `test_results` directory.
