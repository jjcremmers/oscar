import h5py

f = h5py.File('PatchTest4_1.h5', 'r')

print(f.keys())

grp = f['elements']
dset = grp['connectivity']

data = dset[:]

print(data)

dset = grp['pointers']

data = dset[:]

print(data)

grp = f['nodes']
dset = grp['coordinates']

data = dset[:]

print(data)

grp = f['nodeData']
dset = grp['displacements']

data = dset[:]

print(data)

dset = grp['S11']

data = dset[:]

print(data)

print(f.attrs.keys())

x = f.attrs['version']

print(x)

print(f['nodeData'].keys())

